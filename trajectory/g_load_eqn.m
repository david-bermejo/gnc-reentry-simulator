function res = g_load_eqn(AoA, h, p)
    T = p.T_interp(h);
    rho = p.rho_interp(h);
    M = 2000 ./ sqrt(p.gamma .* p.Rg .* T);
    CL = p.lift_clean(AoA, M);
    CD = p.drag_clean(AoA, M);

    %res = sqrt(2.*p.m.*p.g0.*p.n_max ./ (rho .* p.Sref * sqrt(CL.^2 + CD.^2)));
    res = sqrt(2.*p.m.*p.g0.*p.n_max ./ (rho .* p.Sref * sqrt(CL.*cos(deg2rad(AoA)) + CD.*sin(deg2rad(AoA)))));
end