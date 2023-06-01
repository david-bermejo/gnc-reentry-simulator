function [sigma, delta_b] = bank_angle_cmd(x, u_prev, p)
    h = x(1) - p.Rp;
    T = p.T_interp(h);
    rho = p.rho_interp(h);
    g = gravity(x(1));

    M = x(4) / sqrt(p.gamma*p.Rg*T);
    AoA = AoA_curve(M);

    Cm0 = p.pitch_clean(AoA, M);
    delta_b = fzero(@(z) Cm0 + p.pitch_flap(AoA,z,M), u_prev(2));
    if delta_b <= 1e-4
        delta_b = u_prev(2);
    end

    q_inf = 0.5*rho*x(4)^2;
    CD0 = p.drag_clean(AoA, M);
    CL0 = p.lift_clean(AoA, M);
    CDb = p.drag_flap(AoA, delta_b, M);
    CLb = p.lift_flap(AoA, delta_b, M);
    CD = CD0 + CDb;
    CL = CL0 + CLb;
    %L = q_inf*CL*p.Sref;
    
    V_qdot = (p.r*p.Qdot_max/(p.Kq*sqrt(rho/p.Rn)))^(1/3);
    V_n = sqrt(2*p.m*p.g0*p.r*p.n_max/(rho*p.Sref*sqrt(CD^2 + CL^2)));
    V_q = sqrt(2*p.r*p.q_max/rho);
    V_cmd = min([V_qdot, V_n, V_q]);

    %bank_fcn = @(V) 2*p.m/(rho*CL*p.Sref) * (g/V^2 - 1/x(1))*cos(x(5));

    sd = sin(x(3));
    cd = cos(x(3));
    sg = sin(x(5));
    cg = cos(x(5));
    sc = sin(x(6));
    cc = cos(x(6));

    bank_fcn = @(V) 2*p.m/(rho*V^2*CL*p.Sref) * ((g - V^2/x(1))*cg - ...
        2*p.omega_cb*V*cd*sc - p.omega_cb^2*x(1)*cd*(cd*cg + sg*sd*cc));

    %L = 0.5*rho*V_n^2*CL*p.Sref;
    %bank_fcn = @(V) p.m/L * (g - V^2/x(1))*cos(x(5));

    %cb = max([bank_fcn(V_qdot), bank_fcn(V_n), bank_fcn(V_q)]);
    cb = bank_fcn(V_cmd);
    tmp = min(max(cb, -1.0), 1.0);
    
    d_chi = delta_heading(x, p);
    d_chi_max = d_chi_curve(p.d_chi0, p.d_chif, p.M0, p.Mf, M);
    sigma = acos(tmp) * sgn(d_chi, d_chi_max, u_prev(1));
    sigma = max(min(sigma, p.sigma_max), -p.sigma_max);
end