function solution = direct_colocation(problem)
    % Perform direct collocation using trapezoidal rule with surrogate optimization as internal NLP solver
    
    % Inputs:
    % - problem.dynamics: function handle for system dynamics (should take in three arguments, t, x and u)
    % - tspan: time span of the trajectory ([t_initial, t_final])
    % - x0: initial state (column vector)
    % - xN: final state (column vector)
    % - sol_guess: solution guess, improves solver performance (struct with fields t, x, u)
    % - path_bounds: bounds on path variables (matrix with nRows = (nStates + nControls), nCols = 2, with each row in the form [lowerBound, upperBound])
    % - lin_con: linear constraints on the decision variables (struct with fields Aeq and beq for equality constraints, and A and b for inequality constraints)
    % - nonlin_con: nonlinear constraints on the decision variables (function handle, t, x and u)
    
    % Outputs:
    % - t: time vector (column vector)
    % - x: state variables at each time (column vector)
    % - u: control variables at each time (column vector)
    % - fval: objective function value at optimal solution
    % - exitflag: solver exit flag (1 = converged, 0 = failed)

    if numel(unique([...
            length(problem.bounds.x0.low), ...
            length(problem.bounds.x0.upp), ...
            length(problem.bounds.xf.low), ...
            length(problem.bounds.xf.upp), ...
            length(problem.bounds.x.low), ...
            length(problem.bounds.x.upp), ...
            size(problem.guess.x, 1)])) ~= 1
        error('Incorrect number of state variables');
    end

    if numel(unique([...
            length(problem.bounds.u.low), ...
            length(problem.bounds.u.upp), ...
            size(problem.guess.u, 1)])) ~= 1
        error('Incorrect number of control variables');
    end

    nStates = length(problem.bounds.x0.low);
    nControls = length(problem.bounds.u.low);

    if ~isfield(problem.options, 'nIntervals')
        error('Number of discretization intervals not selected.');
    end

    guess.t = problem.guess.t;
    guess.x = problem.guess.x;
    guess.u = problem.guess.u;

    for k=1:length(problem.options)
        nIntervals = problem.options(k).nIntervals; % Number of intervals for collocation
        nDecisionVars = (nStates + nControls) * (nIntervals + 1) + 2;
    
        %% Construction of the initial guess
        ts = linspace(guess.t(1), guess.t(end), nIntervals+1);
    
        xs = zeros(nStates, nIntervals+1);
        for i=1:nStates
            xs(i,:) = interp1(guess.t, guess.x(i,:), ts, "linear");
        end
    
        us = zeros(nControls, nIntervals+1);
        for i=1:nControls
            us(i,:) = interp1(guess.t, guess.u(i,:), ts, "linear");
        end
    
        ys = zeros(nDecisionVars, 1);
        ys(1:end-2) = reshape([xs; us], [], 1);
        ys(end-1) = guess.t(1);
        ys(end) = guess.t(end);
    
        %% Lower and upper bound constraints
        lb = zeros(nDecisionVars, 1);
        ub = zeros(nDecisionVars, 1);
    
        nVars = nStates + nControls;
        lb(1:nVars) = [problem.bounds.x0.low; problem.bounds.u.low];
        lb(nVars+1:end-(nVars+2)) = repmat([problem.bounds.x.low; problem.bounds.u.low], nIntervals-1, 1);
        lb(end-(nVars+1):end-2) = [problem.bounds.xf.low; problem.bounds.u.low];
        lb(end-1) = problem.bounds.t0.low;
        lb(end) = problem.bounds.tf.low;
    
        ub(1:nVars) = [problem.bounds.x0.upp; problem.bounds.u.upp];
        ub(nVars+1:end-(nVars+2)) = repmat([problem.bounds.x.upp; problem.bounds.u.upp], nIntervals-1, 1);
        ub(end-(nVars+1):end-2) = [problem.bounds.xf.upp; problem.bounds.u.upp];
        ub(end-1) = problem.bounds.t0.upp;
        ub(end) = problem.bounds.tf.upp;
        
        %% Solve the NL Problem
        params.nIntervals = nIntervals;
        params.nStates = nStates;
    
        % Set solver and pass options to the internal algorithm
        if isfield(problem.options(k), 'solver')
            solver = problem.options(k).solver;
        else
            solver = 'fmincon';
        end
        
        options = optimoptions( ...
            solver, ...
            'Display', 'iter', ...
            'MaxIterations', 400, ...
            'MaxFunctionEvaluations', 250000);
    
        fn = fieldnames(problem.options(k).opt);
        for i=1:numel(fn)
            options.(fn{i}) = problem.options(k).opt.(fn{i});
        end
    
        % Select quadrature scheme for direct collocation
        if isfield(problem.options(k), 'scheme')
            scheme = problem.options(k).scheme;
        else
            scheme = 'trapezoid';
        end
    
        switch scheme
            case 'trapezoid'
                quad            = @trapz_quad;
                state_interp    = @trapz_state_interp;
                control_interp  = @trapz_control_interp;
                error_eval      = @trapz_error_eval;
            case 'simpson'
                quad            = @simpson_quad;
                state_interp    = @simpson_state_interp;
                control_interp  = @simpson_control_interp;
                error_eval      = @simpson_error_eval;
            otherwise
                error("Invalid quadrature scheme selected. Valid options are: 'trapezoid', 'simpson'.");
        end
    
        [x, fval, exitflag] = fmincon( ...
            @(x) objfun(x, problem.bndObj, problem.pathObj, quad, params), ...
            ys, ...
            [], [], [], [], ...
            lb, ub, ...
            @(x) nlcon(x, problem.pathCst, problem.dynamics, quad, params), ...
            options);

        vars = reshape(x(1:end-2), [], nIntervals+1);
        guess.t = linspace(x(end-1), x(end), nIntervals+1);
        guess.x = vars(1:nStates, :);
        guess.u = vars(nStates+1:end, :);
        guess.fval = fval;
        guess.exitflag = exitflag;
    end

    solution.t = guess.t;
    solution.x = guess.x;
    solution.u = guess.u;
    solution.fval = guess.fval;
    solution.exitflag = guess.exitflag;

    solution.interp.x = @(t) state_interp(t, ...
        solution.t, solution.x, solution.u, problem.dynamics);
    solution.interp.u = @(t) control_interp(t, ...
        solution.t, solution.u);
    solution.interp.err = @(t) error_eval(t, ...
        solution.t, solution.x, solution.u, problem.dynamics);
end