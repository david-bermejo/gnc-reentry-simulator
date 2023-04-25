function out = nav()
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

    sc_ng_acc_eci_init  = [0; 0; 0];
    sc_ang_vel_bod_init = [0; 0; 0];
    b_acc_init          = [0; 0; 0];
    b_gyr_init          = [0; 0; 0];

    out.x0 = [sc_pos_eci_init; sc_vel_eci_init; sc_ng_acc_eci_init; ...
              sc_ang_vel_bod_init; out.sc_att_q_bod2eci_init; ...
              b_acc_init; b_gyr_init];

    % Initial state covariance matrix
    %out.P0 = zeros(22);
    out.Q = zeros(22);
    out.Q(1:3,1:3)     = 1e-1 * eye(3);
    out.Q(4:6,4:6)     = 1e-1 * eye(3);
    out.Q(7:9,7:9)     = 1e-6 * eye(3);
    out.Q(10:12,10:12) = 1e-9 * eye(3);
    out.Q(13:16,13:16) = 1e-9 * eye(4);
    out.Q(17:19,17:19) = 1e-6 * eye(3);
    out.Q(20:22,20:22) = 1e-6 * eye(3);
    out.P0 = out.Q;

    out.R = zeros(6);
    out.R(1:3,1:3) = 1e-8 * eye(3);
    out.R(4:6,4:6) = 1e-8 * eye(3);

    out.dq_K = 1.5;
    out.Ma = zeros(3);
    out.Mo = zeros(3);
end