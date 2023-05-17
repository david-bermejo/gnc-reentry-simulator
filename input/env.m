function out = env()
    %% Constants
    % Planet radius:
    out.Rp       = 3396.2 * 1e3;    % [m]
    % Planet angular velocity:
    out.omega_cb = 7.07765809e-5;   % [rad/s]
    % Planet oblateness:
    out.f        = 0.00589;         % [-]
    % Planet gravitational constant:
    out.mu       = 0.042828e15;     % [m^3/s^2]
    % Standard Gravity:
    out.g0       = 9.80665;         % [m/s^2]
    % Gas constant:
    out.gamma    = 1.2941;          % [-]
    % Gas constant:
    out.Rg       = 8314.3/44.01;    % [-]

    %% Atmospheric Model (Mars Climate Database)
    out = append2struct(out, {'h_base', 'T_base', 'rho_base'}, atmospheric_coeffs(), 3);

    %% Initial Conditions (DKE propagator)
    load('trajectory.mat', 'tau0', 'delta0', 'V0', 'gamma0', 'chi0');
    % Initial radial position:
    R0      = out.Rp + 120e3;   % [m]
    % Initial angle of attack:
    alpha0  = deg2rad(40);      % [deg]
    % Initial angle of sideslip:
    beta0   = 0;                % [deg]
    % Initial commanded bank angle:
    sigma0  = deg2rad(-165);    % [deg]
    % Initial attitude quaternion:
    q0      = dcm2quat(...
                z2dcm(-tau0) * y2dcm(pi/2 + delta0) * ...
                z2dcm(-chi0) * y2dcm(-gamma0) * ...
                x2dcm(sigma0) * z2dcm(beta0) * y2dcm(-alpha0));

    % Initial state (16x1):
    out.x0  = [R0; tau0; delta0; ...
               V0; gamma0; chi0; ...
               0; 0; 0; ...
               alpha0; beta0; sigma0; ...
               q0];
end