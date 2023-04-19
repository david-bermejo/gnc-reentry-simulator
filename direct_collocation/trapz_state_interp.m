function x = trapz_state_interp(t, ts, xs, us, dynamics)
    if t(1) < ts(1) || t(end) > ts(end)
        error("Time interpolation out of bounds.");
    end

    x = zeros(size(xs, 1), length(t));
    fs = dynamics(ts, xs, us);
    idx = 1;

    for i=1:length(t)
        while t(i) > ts(idx+1)
            idx = idx + 1;
        end

        tau = t(i) - ts(idx);
        h = ts(idx+1) - ts(idx);
        x(:,i) = xs(:,idx) + fs(:,idx).*tau + tau.^2./(2.*h) .* (fs(:,idx+1) - fs(:,idx));
    end
end