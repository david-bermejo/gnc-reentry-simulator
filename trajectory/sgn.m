function res = sgn(d_chi, chi_max, u)
    res = ones(size(u));

    if length(u) == 1
        if d_chi >= chi_max
            res = -1;
        elseif d_chi <= -chi_max
            res = 1;
        else
            res = u/(abs(u) + (u == 0));
        end
    else
        for i=2:length(u)
            if d_chi(i) >= chi_max
                res(i) = -1;
            elseif d_chi(i) <= -chi_max
                res(i) = 1;
            else
                res(i) = res(i-1);
            end
        end
    end
end