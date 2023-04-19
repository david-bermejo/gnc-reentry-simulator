function res = trapz_quad(ts, xs, us, fcn)
    hk = ts(2:end) - ts(1:end-1);
    fs = fcn(ts, xs, us);
    res = hk./2 .* (fs(:,1:end-1) + fs(:,2:end));
end