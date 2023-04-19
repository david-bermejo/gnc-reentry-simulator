function res = simpson_quad(ts, xs, us, fcn)
    hk = ts(2:end) - ts(1:end-1);
    fs = fcn(ts, xs, us);
    tmid = (ts(2:end) + ts(1:end-1))./2;
    xmid = (xs(:,1:end-1) + xs(:,2:end))./2 + ...
        (ts(2:end) - ts(1:end-1))./8 .* (fs(:,1:end-1) - fs(:,2:end));
    umid = (us(:,1:end-1) + us(:,2:end)) ./ 2;
    fmid = fcn(tmid, xmid, umid);
    
    res = hk./6 .* (fs(:,1:end-1) + 4.*fmid + fs(:,2:end));
end