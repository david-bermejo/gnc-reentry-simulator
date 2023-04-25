function out = env()
    out.omega_cb = 7.07765809e-5;     % [rad/s]
    out.mu = 0.042828e15;             % [m^3/s^2]
    out.Rp = 3396.2 * 1e3;            % [m]
    out.f = 0.00589;                  % [-]
    out.gamma = 1.2941;               % [-]
    out.Rg = 8314.3/44.01;            % [-]
    out.g0 = 9.80665;                 % [m/s^2]

    % Atmosphere model (Mars Climate Database)
    out = append2struct(out, {'h_base', 'T_base', 'rho_base'}, atmospheric_coeffs(), 3);

    % Initial Conditions
    load('trajectory.mat', 'tau0', 'delta0', 'V0', 'gamma0', 'chi0');
    R0      = out.Rp + 120e3;
    alpha0  = deg2rad(40);
    beta0   = 0;
    sigma0  = deg2rad(-165);
    q0      = dcm2quat(...
                z2dcm(-tau0) * y2dcm(pi/2 + delta0) * ...
                z2dcm(-chi0) * y2dcm(-gamma0) * ...
                x2dcm(sigma0) * z2dcm(beta0) * y2dcm(-alpha0));

    out.x0  = [R0; tau0; delta0; ...
               V0; gamma0; chi0; ...
               0; 0; 0; ...
               alpha0; beta0; sigma0; ...
               q0];
    out.e0 = V0^2/2 + out.mu/R0^2 * (R0 - out.Rp);
end