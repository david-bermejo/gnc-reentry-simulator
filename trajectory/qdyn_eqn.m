function res = qdyn_eqn(h, p)
    rho = p.rho_interp(h);
    res = sqrt(2.*p.q_max ./ rho);
end