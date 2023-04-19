function u = simpson_control_interp(t, ts, us)
    if t(1) < ts(1) || t(end) > ts(end)
        error("Time interpolation out of bounds.");
    end

    u = zeros(size(us, 1), length(t));
    umid = (us(:,1:end-1) + us(:,2:end)) ./ 2;

    idx = 1;

    for i=1:length(t)
        while t(i) > ts(idx+1)
            idx = idx + 1;
        end

        tau = t(i) - ts(idx);
        h = ts(idx+1) - ts(idx);

        u(:,i) = 2./(h.^2).*(tau-h./2).*(tau-h).*us(:,idx) - ...
            4./(h.^2).*tau.*(tau-h).*umid(:,idx) + ...
            2./(h.^2).*tau.*(tau-h./2).*us(:,idx+1);
    end
end