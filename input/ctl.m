function out = ctl(noiseFlag)
    %% Time definitions
    % Execution frequency:
    out.freq = 20;          % [Hz]
    % Execution time interval:
    out.tsamp = 1/out.freq; % [s]

    %% Longitudinal LQR Controller configuration
    % Maximum pitch rate deviation:
    q_max       = deg2rad(50);  % [deg/s]
    % Maximum angle of attack deviation:
    AoA_max     = deg2rad(0.5); % [deg]
    % Maximum commanded elevator deflection:
    delta_e_max = 40;       % [deg]
    % Maximum commanded RCS torque (y axis):
    Ty_max      = 10400;    % [Nm]

    % State cost matrix:
    out.Q_lon      = zeros(2);
    out.Q_lon(1,1) = 1/q_max^2;
    out.Q_lon(2,2) = 1/AoA_max^2;

    % Control cost matrix:
    out.R_lon      = zeros(2);
    out.R_lon(1,1) = 1/delta_e_max^2;
    out.R_lon(2,2) = 1/Ty_max^2;

    %% Lateral LQR Controller configuration
    % Maximum roll rate deviation:
    p_max       = deg2rad(50);  % [deg/s]
    % Maximum yaw rate deviation:
    r_max       = deg2rad(50);  % [deg/s]
    % Maximum angle of sideslip deviation:
    AoS_max     = deg2rad(2); % [deg]
    % Maximum bank angle deviation:
    sigma_max   = deg2rad(8);   % [deg]
    % Maximum commanded aileron deflection:
    delta_a_max = deg2rad(40);  % [deg]
    % Maximum commanded rudder deflection:
    delta_r_max = deg2rad(40);  % [deg]
    % Maximum commanded RCS torque (x axis):
    Tx_max      = 1600;         % [Nm]
    % Maximum commanded RCS torque (z axis):
    Tz_max      = 7600;         % [Nm]

    % State cost matrix:
    out.Q_lat      = zeros(4);
    out.Q_lat(1,1) = 1/p_max^2;
    out.Q_lat(2,2) = 1/r_max^2;
    out.Q_lat(3,3) = 1/AoS_max^2;
    out.Q_lat(4,4) = 1/sigma_max^2;

    % Control cost matrix:
    out.R_lat      = zeros(4);
    out.R_lat(1,1) = 1/delta_a_max^2;
    out.R_lat(2,2) = 1/delta_r_max^2;
    out.R_lat(3,3) = 1/Tx_max^2;
    out.R_lat(4,4) = 1/Tz_max^2;
end