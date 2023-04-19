function eps = trapz_error_eval(t, ts, xs, us, dynamics)
    if t(1) < ts(1) || t(end) > ts(end)
        error("Time interpolation out of bounds.");
    end

    eps = zeros(size(xs, 1), length(t));
    fs = dynamics(ts, xs, us);
    idx = 1;

    for i=1:length(t)
        while t(i) > ts(idx+1)
            idx = idx + 1;
        end

        tau = t(i) - ts(idx);
        h = ts(idx+1) - ts(idx);

        xdot = fs(:,idx) + (fs(:,idx+1) - fs(:,idx)) .* (tau./h);
        x = xs(:,idx) + fs(:,idx).*tau + tau.^2./(2.*h) .* (fs(:,idx+1) - fs(:,idx));
        u = us(:,idx) + (us(:,idx+1) - us(:,idx)) .* (tau./h);

        f = dynamics(t(i), x, u);
        eps(:,i) = xdot - f;
    end
end