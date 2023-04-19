function [c, ceq] = nlcon(x, pathCst, dynamics, quadrature, p)
    ts = linspace(x(end-1), x(end), p.nIntervals+1);
    vars = reshape(x(1:end-2), [], p.nIntervals+1);
    states = vars(1:p.nStates, :);
    controls = vars(p.nStates+1:end, :);

    [Cin, Ceq] = pathCst(ts, states, controls);
    
    dyn_eq = reshape(states(:,2:end) - states(:,1:end-1) - ...
        quadrature(ts, states, controls, dynamics), [], 1);

    c = [reshape(Cin, [], 1)];
    ceq = [reshape(Ceq, [], 1); dyn_eq];
end