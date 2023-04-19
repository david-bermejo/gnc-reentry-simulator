function [value, isterminal, direction] = myEvent(t, x, p)
    value      = x(1) - p.xf(1);
    isterminal = 1;   % Stop the integration
    direction  = -1;
end