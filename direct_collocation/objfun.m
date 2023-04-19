function Fval = objfun(x, bndObj, pathObj, quadrature, p)
    ts = linspace(x(end-1), x(end), p.nIntervals+1);
    vars = reshape(x(1:end-2), [], p.nIntervals+1);
    states = vars(1:p.nStates, :);
    controls = vars(p.nStates+1:end, :);

    Fval = bndObj(states(:,1), ts(1), states(:,end), ts(end)) + ...
        sum(quadrature(ts, states, controls, pathObj));
end