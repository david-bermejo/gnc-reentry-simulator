function out = sc()
    % Geometry definition constants
    out.Sref = 110;                     % [m^2]
    out.bref = 13;                      % [m]
    out.cref = 23;                      % [m]
    out.m = 26029.0;                    % [kg]
    
    % Spacecraft inertia matrix in BRF [kg/s^2]:
    out.I = [119605,      0, -20372;
                  0, 769000,      0;
             -20372,      0, 805395];

    % RCS Thrusters and Aero Surfaces data:
    out.RCS_T_on = 400;                 % [N]
    out.delta_r_bounds = [0, 40];       % [deg]
    out.delta_e_bounds = [-40, 40];     % [deg]
    out.delta_b_bounds = [-20, 30];     % [deg]
    out.x_cp = [0; 0; 0];               % [m]

    out.x_com = [13; 0; 0];
    out.x_RCS_roll_forward = [20; 2; 0];    % [m]
    out.x_RCS_roll_rear = [20; -2; 0];      % [m]
    out.x_RCS_pitch_forward = [1; 0; 0];    % [m]
    out.x_RCS_pitch_rear = [20; 0; 0];      % [m]
    out.x_RCS_yaw_forward = [1; 0; 0.3];    % [m]
    out.x_RCS_yaw_rear = [20; 0; 0];        % [m]
end
