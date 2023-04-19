function res = delta_heading(x, p)
    X = cos(p.xf(3)) .* sin(p.xf(2) - x(2,:));
    Y = cos(x(3,:)) .* sin(p.xf(3)) - sin(x(3,:)) .* cos(p.xf(3)) .* cos(p.xf(2) - x(2,:));
    chi_des = atan2(X, Y);
    res = chi_des - x(6,:);
end