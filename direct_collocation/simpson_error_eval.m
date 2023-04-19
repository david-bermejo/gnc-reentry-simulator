function eps = simpson_error_eval(t, ts, xs, us, dynamics)
    if t(1) < ts(1) || t(end) > ts(end)
        error("Time interpolation out of bounds.");
    end

    eps = zeros(size(xs, 1), length(t));
    fs = dynamics(ts, xs, us);

%     xmid = (xs(:,1:end-1) + xs(:,2:end))./2 + ...
%         (ts(2:end) - ts(1:end-1))./8 .* (fs(:,1:end-1) - fs(:,2:end));
    fmid = -3./(2 .* (ts(2:end) - ts(1:end-1))) .* (xs(:,1:end-1) - xs(:,2:end)) - ...
        (fs(:,1:end-1) + fs(:,2:end))./4;
    umid = (us(:,1:end-1) + us(:,2:end)) ./ 2;

    idx = 1;

    for i=1:length(t)
        while t(i) > ts(idx+1)
            idx = idx + 1;
        end

        tau = t(i) - ts(idx);
        h = ts(idx+1) - ts(idx);
        x = xs(:,idx) + fs(:,idx) .* tau + ...
            (-3.*fs(:,idx) + 4.*fmid(:,idx) - fs(:,idx+1)) .* tau.^2./(2.*h) + ...
            (2.*fs(:,idx) - 4.*fmid(:,idx) + 2.*fs(:,idx+1)) .* tau.^3./(3.*h.^2);
        
        u = 2./(h.^2).*(tau-h./2).*(tau-h).*us(:,idx) - ...
            4./(h.^2).*tau.*(tau-h).*umid(:,idx) + ...
            2./(h.^2).*tau.*(tau-h./2).*us(:,idx+1);

        xdot = 2./(h.^2).*(tau-h./2).*(tau-h).*fs(:,idx) - ...
            4./(h.^2).*tau.*(tau-h).*fmid(:,idx) + ...
            2./(h.^2).*tau.*(tau-h./2).*fs(:,idx+1);

        f = dynamics(t(i), x, u);
        eps(:,i) = xdot - f;
    end
end