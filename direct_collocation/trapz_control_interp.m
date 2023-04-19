function u = trapz_control_interp(t, ts, us)
    if t(1) < ts(1) || t(end) > ts(end)
        error("Time interpolation out of bounds.");
    end

    u = zeros(size(us, 1), length(t));
    idx = 1;

    for i=1:length(t)
        while t(i) > ts(idx+1)
            idx = idx + 1;
        end

        tau = t(i) - ts(idx);
        h = ts(idx+1) - ts(idx);
        u(:,i) = us(:,idx) + (us(:,idx+1) - us(:,idx)) .* (tau./h);
    end
end