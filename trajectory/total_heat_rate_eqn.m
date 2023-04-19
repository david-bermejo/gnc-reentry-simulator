function res = total_heat_rate_eqn(V, h, p)
    rho = p.rho_interp(h);
    q_conv = 1.9027e-4 .* sqrt(rho./p.Rn) .* V.^3;
    res = q_conv ./ 1000 - p.Qdot_max;
end