function out = nav_w_bias()
    out.freq = 20;
    out.tsamp = 1/out.freq;

    % Initial Conditions
    load('trajectory.mat', 'tau0', 'delta0', 'V0', 'gamma0', 'chi0');
    R0      = (3396.2 + 120) * 1e3;
    alpha0  = deg2rad(40);
    beta0   = 0;
    sigma0  = deg2rad(-165);
    
    % Initial state
    out.sc_pos_vrt_init         = [R0; tau0; delta0];
    out.sc_vel_tg_init          = [V0; gamma0; chi0];
    out.sc_ang_vel_bod_init     = [0; 0; 0];
    out.sc_ae_ang_gnd_init      = [alpha0; beta0; sigma0];
    out.sc_att_q_bod2eci_init   = dcm2quat(...
                                  z2dcm(-tau0) * y2dcm(pi/2 + delta0) * ...
                                  z2dcm(-chi0) * y2dcm(-gamma0) * ...
                                  x2dcm(sigma0) * z2dcm(beta0) * y2dcm(-alpha0));

    % Initial state
    sc_pos_eci_init     = [R0*cos(tau0)*cos(delta0);
                           R0*sin(tau0)*cos(delta0);
                                     R0*sin(delta0)];

    sc_vel_vrt_init     = [V0*cos(chi0)*cos(gamma0);
                           V0*sin(chi0)*cos(gamma0);
                                    -V0*sin(gamma0)];

    omega_cb_vec = [0; 0; 7.07765809e-5]; % [rad/s]
    sc_vel_eci_init     = z2dcm(-tau0)*y2dcm(pi/2 + delta0) * ...
                          sc_vel_vrt_init + cross(omega_cb_vec, sc_pos_eci_init);

    out.x0 = [sc_pos_eci_init; sc_vel_eci_init; ...
              zeros(3,1); out.sc_att_q_bod2eci_init; zeros(3,1)];

    %% Extended Kalman Filter Configuration
    out.x_nstates       = 16;
    out.x_sta_pos_idx   = 1:3;
    out.x_sta_vel_idx   = 4:6;
    out.x_sta_ba_idx    = 7:9;
    out.x_sta_att_idx   = 10:13;
    out.x_sta_bw_idx    = 14:16;

    out.dx_nstates      = 15;
    out.dx_sta_pos_idx  = 1:3;
    out.dx_sta_vel_idx  = 4:6;
    out.dx_sta_ba_idx   = 7:9;
    out.dx_sta_att_idx  = 10:12;
    out.dx_sta_bw_idx   = 13:15;

    %% Neccesary data
    % Gravity Reference
    g0 = 9.80665; % [m/s^2]

    % Acceleration Bias [m/s^2]
    acc_bias_var    = (0.1*1e-6*g0)^2; % [ug]

    % Velocity Random Walk [m/s^(3/2)]
    N_VRW           = 0.1; % [ug/sqrt(Hz)]
    vrw_var         = (N_VRW*1e-6*g0)^2;

    % Acceleration Random Walk [m/s^(5/2)]
    N_ARW           = 0.05; % [ug/sqrt(hr)]
    xrw_var         = (N_ARW*1e-6*g0/60)^2;

    % Gyro Bias [rad/s]
    gyr_bias_var    = (deg2rad(0.01)/3600)^2; % [deg/h]

    % Angular Random Walk [rad/sqrt(s)]
    N_ARW           = 0.005; % [deg/sqrt(hr)]
    arw_var         = (deg2rad(N_ARW)/60)^2;

    % Rate Random Walk [rad/s^2]
    N_RRW           = 0.005; % [deg/hr^(3/2)]
    rrw_var         = (deg2rad(N_RRW)/60^3)^2;

    %% Initial state covariance matrix
    out.Q = zeros(12);
    out.Q(1:3,1:3)     = vrw_var * eye(3);
    out.Q(4:6,4:6)     = xrw_var * eye(3);
    out.Q(7:9,7:9)     = arw_var * eye(3);
    out.Q(10:12,10:12) = rrw_var * eye(3);

    out.P0 = zeros(15);
    out.P0(1:3,1:3)     = 100^2 * eye(3);
    out.P0(4:6,4:6)     = 0.1^2 * eye(3);
    out.P0(7:9,7:9)     = acc_bias_var * eye(3);
    out.P0(10:12,10:12) = (1e-4)^2 * eye(3);
    out.P0(13:15,13:15) = gyr_bias_var * eye(3);
    out.P0 = out.P0 * 1000;

    out.R = zeros(10);
    out.R(1:5,1:5) = 1e1 * eye(5);
    out.R(6:10,6:10) = 1e-2 * eye(5);

    out.dq_K = 1.5;
end