function out = nav()
    %% Constants
    % Planet radius:
    Rp = 3396.2 * 1e3;          % [m]
    % Planet angular velocity:
    omega_cb = 7.07765809e-5;   % [rad/s]
    % Planet rotation vector:
    omega_cb_vec = [0; 0; omega_cb];
    % Standard Gravity:
    g0 = 9.80665;               % [m/s^2]

    %% Time definitions
    % Sample frequency:
    out.freq = 20;          % [Hz]
    % Sample time:
    out.tsamp = 1/out.freq; % [s]

    %% SC initial conditions (offline trajectory)
    load('trajectory-v33.mat', 'tau0', 'delta0', 'V0', 'gamma0', 'chi0', 'u');
    % Initial radial position:
    R0      = Rp + 120e3;       % [m]
    % Initial angle of attack:
    alpha0  = deg2rad(45);      % [deg]
    % Initial angle of sideslip:
    beta0   = 0;                % [deg]
    % Initial bank angle:
    sigma0  = u(1,1);           % [deg]
    
    % Initial radial position, longitude & latitude:
    out.sc_pos_vrt_init         = [R0; tau0; delta0];
    % Initial velocity, glideslope & heading:
    out.sc_vel_tg_init          = [V0; gamma0; chi0];
    % Initial angular velocity (body w.r.t. inertial):
    out.sc_ang_vel_bod_init     = [0; 0; 0];
    % Initial aerodynamic angles:
    out.sc_ae_ang_gnd_init      = [alpha0; beta0; sigma0];
    % Initial attitude quaternion (body w.r.t. inertial):
    out.sc_att_q_bod2eci_init   = dcm2quat(...
                                  z2dcm(-tau0) * y2dcm(pi/2 + delta0) * ...
                                  z2dcm(-chi0) * y2dcm(-gamma0) * ...
                                  x2dcm(sigma0) * z2dcm(beta0) * ...
                                  y2dcm(-alpha0));

    %% Extended Kalman Filter initial state
    % Initial position (inertial frame):
    sc_pos_eci_init     = [R0*cos(tau0)*cos(delta0);
                           R0*sin(tau0)*cos(delta0);
                                     R0*sin(delta0)];
    % Initial velocity (vertical frame):
    sc_vel_vrt_init     = [V0*cos(chi0)*cos(gamma0);
                           V0*sin(chi0)*cos(gamma0);
                                    -V0*sin(gamma0)];
    % Initial velocity (inertial frame):
    sc_vel_eci_init     = z2dcm(-tau0)*y2dcm(pi/2 + delta0) * ...
                          sc_vel_vrt_init + cross(omega_cb_vec, sc_pos_eci_init);

    % Initial EKF state (pos, vel, att_q):
    out.x0 = [sc_pos_eci_init; sc_vel_eci_init; out.sc_att_q_bod2eci_init];

    %% Extended Kalman Filter Configuration
    % State indices
    out.x_nstates       = 10;
    out.x_sta_pos_idx   = 1:3;
    out.x_sta_vel_idx   = 4:6;
    out.x_sta_att_idx   = 7:10;

    % State increment indices
    out.dx_nstates      = 9;
    out.dx_sta_pos_idx  = 1:3;
    out.dx_sta_vel_idx  = 4:6;
    out.dx_sta_att_idx  = 7:9;

    % Initial state error covariance matrix
    out.P0 = zeros(9);
    out.P0(1:3,1:3) = 500^2 * eye(3);
    out.P0(4:6,4:6) = 10^2 * eye(3);
    out.P0(7:9,7:9) = (1e-3)^2 * eye(3);
    out.P0 = out.P0;
    
    % Velocity Random Walk [m/s^(3/2)]
    N_VRW      = 0.1;   % [mg/sqrt(Hz)]
    vrw_var    = (N_VRW*1e-3*g0)^2;

    % Angular Random Walk [rad/sqrt(s)]
    N_ARW      = 0.005; % [deg/sqrt(hr)]
    arw_var    = (deg2rad(N_ARW)/60)^2;

    % Process noise covariance matrix
    out.Q = zeros(6);
    out.Q(1:3,1:3)  = vrw_var * eye(3);
    out.Q(4:6,4:6)  = arw_var * eye(3);

    % Process noise covariance matrix (alt.)
    out.Q2 = zeros(9);
    out.Q2(1:3,1:3) = 1^2 * eye(3);
    out.Q2(4:6,4:6) = 0.01^2 * eye(3);
    out.Q2(7:9,7:9) = (1e-5)^2 * eye(3);
    out.Q2 = out.Q2;

    % Measurement noise covariance matrix
    out.R = zeros(8);
    out.R(1:4,1:4) = 10^2 * eye(4);
    out.R(5:8,5:8) = 0.1^2 * eye(4);
end