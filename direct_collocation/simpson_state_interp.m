function x = simpson_state_interp(t, ts, xs, us, dynamics)
    if t(1) < ts(1) || t(end) > ts(end)
        error("Time interpolation out of bounds.");
    end

    x = zeros(size(xs, 1), length(t));
    fs = dynamics(ts, xs, us);

    fmid = -3./(2 .* (ts(2:end) - ts(1:end-1))) .* (xs(:,1:end-1) - xs(:,2:end)) - ...
        (fs(:,1:end-1) + fs(:,2:end))./4;

    idx = 1;

    for i=1:length(t)
        while t(i) > ts(idx+1)
            idx = idx + 1;
        end

        tau = t(i) - ts(idx);
        h = ts(idx+1) - ts(idx);

        x(:,i) = xs(:,idx) + fs(:,idx) .* tau + ...
            (-3.*fs(:,idx) + 4.*fmid(:,idx) - fs(:,idx+1)) .* tau.^2./(2.*h) + ...
            (2.*fs(:,idx) - 4.*fmid(:,idx) + 2.*fs(:,idx+1)) .* tau.^3./(3.*h.^2);
    end
end