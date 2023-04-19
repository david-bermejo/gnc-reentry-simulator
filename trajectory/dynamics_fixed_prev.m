function xdot = dynamics_fixed(t, x, u, p)
    % x = [R, lon, lat, V, gamma, chi]
    % u = [sigma]
    
    sd = sin(x(3,:));
    cd = cos(x(3,:));
    td = tan(x(3,:));
    sg = sin(x(5,:));
    cg = cos(x(5,:));
    sc = sin(x(6,:));
    cc = cos(x(6,:));
    
    h = x(1,:) - p.Rp;
    T = p.T_interp(h);
    rho = p.rho_interp(h);
    g = gravity(x(1,:));

    M = x(4,:) ./ sqrt(p.gamma.*p.Rg.*T);
    %AoA = ones(size(M)) * 40.0;
    AoA = AoA_curve(M);

    q_inf = 0.5.*rho.*x(4,:).^2;
    CD = p.drag_clean(AoA, M);
    CL = p.lift_clean(AoA, M);
    D = q_inf.*CD.*p.Sref;
    L = q_inf.*CL.*p.Sref;

    Fv = -D - p.m.*g.*sg;
    Fg = L.*cos(u) - p.m.*g.*cg;
    Fc = -L.*sin(u);
    
    xdot = zeros(size(x));
    xdot(1,:) = x(4,:).*sg;
    xdot(2,:) = x(4,:).*sc.*cg ./ (x(1,:).*cd);
    xdot(3,:) = x(4,:).*cc.*cg ./ x(1,:);
    xdot(4,:) = Fv./p.m + p.omega_cb.^2.*x(1,:).*cd.*(sg.*cd - cg.*sd.*cc);
    xdot(5,:) = (Fg./p.m + 2.*p.omega_cb.*x(4,:).*cd.*sc + x(4,:).^2./x(1,:).*cg + ...
        p.omega_cb.^2.*x(1,:).*cd.*(cd.*cg + sg.*sd.*cc)) ./ x(4,:);
    xdot(6,:) = (Fc./p.m + 2.*p.omega_cb.*x(4,:).*(sd.*cg - cd.*sg.*cc) + ...
        x(4,:).^2./x(1,:).*cg.^2.*td.*sc + p.omega_cb.^2.*x(1,:).*cd.*sd.*sc)./(x(4,:).*cg);
end