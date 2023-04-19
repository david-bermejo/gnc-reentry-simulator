function res = glide_eqn(AoA, V, h, p)
    R = p.Rp + h;
    g = gravity(R);
    W = p.m .* g;

    T = p.T_interp(h);
    rho = p.rho_interp(h);
    M = V ./ sqrt(p.gamma .* p.Rg .* T);
    CL = p.lift_clean(AoA, M);

    res = 2.*(W./p.Sref)./CL .* (1./V.^2 - 1./(R.*g)) - rho;
end