function res = dynamics_example()
    addpath(pwd + "\..\fcn\aerodynamic_coefficients")
    addpath(pwd + "\..\fcn\atmospheric_model")
    addpath(pwd + "\..\fcn\gravity_model")

    %% Parameters setup
    params.mu = 0.042828e15;            % [m^3/s^2]
    params.Rp = 3396.2 * 1e3;           % [m]
    params.omega_cb = 7.07765809e-5;    % [rad/s]
    params.gamma = 1.2941;              % [-]
    params.Rg = 8314.3/44.01;           % [-]
    params.Sref = 110;                  % [m^2]
    params.Rn = 1.2;                    % [m]
    params.m = 26029.0;                 % [kg]
    params.g0 = 9.81;                   % [m/s^2]
    
    params.Qdot_max = 450;              % [kW/m^2]
    params.q_max = 5000;                % [Pa]
    params.n_max = 2.5;                 % [-]

    params.drag_clean = drag_clean_interpolant('spline');
    params.lift_clean = lift_clean_interpolant('spline');
    [params.T_interp, params.p_interp, params.rho_interp] = atmosphere('spline');

    %% Initial conditions
    R0 = params.Rp + 120.0e3;           % [m]
    tau0 = deg2rad(10.3);               % [rad]
    delta0 = deg2rad(0);             % [rad]
    V0 = 7000;                          % [m/s]
    gamma0 = deg2rad(-10.5);             % [rad]
    chi0 = deg2rad(33);               % [rad]
    params.x0 = [R0, tau0, delta0, V0, gamma0, chi0]';
    
    %% Terminal conditions
    Rf = params.Rp + 10e3;               % [m]
    tauf = deg2rad(21);               % [rad]
    deltaf = deg2rad(20);             % [rad]
    Vf = 800.0;                         % [m/s]
    gammaf = deg2rad(-2.0);             % [rad]
    chif = deg2rad(90.0);               % [rad]
    params.xf = [Rf, tauf, deltaf, Vf, gammaf, chif]';

    %% Calculate the initial guess
    N = 1000;
    tspan = [0, 500];
    ts = linspace(tspan(1), tspan(2), N);
    options = odeset('MaxStep',2, 'Events', @(t, x) myEvent(t, x, params));
    u = zeros(1,N);
    y = zeros(6,N);
    y(:,1) = params.x0;
    
    for i=2:N
        h = y(1,i-1) - params.Rp;
    
        if h <= 60e3
            T = params.T_interp(h);
            rho = params.rho_interp(h);
            g = gravity(y(1,i-1));
        
            M = y(4,i-1) ./ sqrt(params.gamma.*params.Rg.*T);
            AoA = AoA_curve(M);
        
            q_inf = 0.5.*rho.*y(4,i-1).^2;
            CL = params.lift_clean(AoA, M);
            L = q_inf.*CL.*params.Sref;
        
            var = (g - y(4,i-1).^2./y(1,i-1)) .* cos(y(5,i-1)) .* params.m ./ L;
            tmp = min(max(var, -1.0), 1.0);
            
            d_chi = delta_heading(y(:,i-1), params);
            u(1,i) = acos(tmp) .* sgn(d_chi, deg2rad(5), u(1,i-1));
        end
    
        % Solve the problem
        solution = ode45(@(t, x) dynamics_fixed(t, x, u(1,i), params), [ts(i-1), ts(i)], y(:,i-1), options);
    
        % Save solution
        y(:,i) = solution.y(:,end);
    
        if y(1,i) <= params.xf(1)
            ts = ts(1:i);
            y = y(:,1:i);
            u = u(1:i);
            break;
        end
    end

    %% Call the function
    profile on
    for i=1:10000
        res = dynamics(ts, y, u, params);
    end
    profile viewer
end