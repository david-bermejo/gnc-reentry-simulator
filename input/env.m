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
    % Gas constant #2:
    out.Rg       = 8314.3/44.01;    % [-]

    %% Atmospheric Model (Mars Climate Database)
    out = append2struct(out, {'h_base', 'T_base', 'rho_base'}, atmospheric_coeffs(), 3);

    % Initial atmospheric density bias (relative percentage):
    out.rho_bias = min(max(randn(1,1)*0.1, -0.2), 0.3);  % [-]
    % Atmospheric density noise variance:
    out.rho_var  = 0.025^2;           % [-]
    % Atmospheric density noise seed:
    out.rho_seed = randi(100000,1); % [-]
    % Atmospheric density autocorrelation time:
    out.rho_Tc   = 100;             % [s]

    %% Initial Conditions (DKE propagator)
    load('trajectory-v33.mat', 'tau0', 'delta0', 'V0', 'gamma0', 'chi0', 'u');
    % Initial radial position:
    R0_mu       = out.Rp + 120e3;   % [m]
    R0_std      = 100;              % [m]
    R0          = randn(1)*R0_std + R0_mu;

    % Initial longitude:
    tau0_mu     = tau0;             % [rad]
    tau0_std    = deg2rad(0.02);    % [rad]
    tau0        = randn(1)*tau0_std + tau0_mu;

    % Initial latitude:
    delta0_mu   = delta0;           % [rad]
    delta0_std  = deg2rad(0.02);    % [rad]
    delta0      = randn(1)*delta0_std + delta0_mu;

    % Initial velocity:
    V0_mu       = V0;               % [m/s]
    V0_std      = 5;               % [m/s]
    V0          = randn(1)*V0_std + V0_mu;
    
    % Initial glideslope:
    gamma0_mu   = gamma0;           % [rad]
    gamma0_std  = deg2rad(0.02);     % [rad]
    gamma0      = randn(1)*gamma0_std + gamma0_mu;

    % Initial heading:
    chi0_mu     = chi0;             % [rad]
    chi0_std    = deg2rad(0.02);     % [rad]
    chi0        = randn(1)*chi0_std + chi0_mu;

    % Initial angle of attack:
    alpha0_mu   = deg2rad(45);      % [rad]
    alpha0_std  = deg2rad(0.1);     % [rad]
    alpha0      = randn(1)*alpha0_std + alpha0_mu;

    % Initial angle of sideslip:
    beta0_mu    = 0;                % [rad]
    beta0_std   = deg2rad(0.1);     % [rad]
    beta0       = randn(1)*beta0_std + beta0_mu;

    % Initial commanded bank angle:
    sigma0_mu   = u(1,1);           % [rad]
    sigma0_std  = deg2rad(0.1);     % [rad]
    sigma0      = randn(1)*sigma0_std + sigma0_mu;

    % Initial attitude quaternion:
    q0          = dcm2quat(...
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