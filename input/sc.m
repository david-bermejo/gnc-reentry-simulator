function out = sc()
    %% Geometry related constants
    % Reference area:
    out.Sref    = 110;      % [m^2]
    % Reference lateral length:
    out.bref    = 13;       % [m]
    % Reference longitudinal length:
    out.cref    = 23;       % [m]
    % SC mass (landing configuration):
    out.m       = 26029.0;  % [kg]
    % SC inertia matrix in Body frame [kg*s^2]:
    out.I = [119605,      0, -20372;
                  0, 769000,      0;
             -20372,      0, 805395];

    %% Aero surfaces data
    % Rudder maximum deflections:
    out.delta_r_bounds  = [0, 40];      % [deg]
    % Elevon maximum deflections:
    out.delta_e_bounds  = [-40, 40];    % [deg]
    % Body flap maximum deflections:
    out.delta_b_bounds  = [-20, 30];    % [deg]
    % Distance between CoP and CoM:
    out.x_cp            = [0; 0; 0];    % [m]

    %% RCS Thrusters data
    % Maximum force exerted by a single thruster:
    out.Fmax    = 400;      % [N]
    % Propellant specific impulse (Hydrazine):
    out.Isp     = 220;      % [s]
    % Standard gravity:
    out.g0      = 9.80665;  % [m/s^2]

    % Position of each thruster (w.r.t. CoM):
    out.RCS_pos = [ % Roll Thrusters (Tx)
                     7.0,  2.0, 0.0;
                     7.0,  2.0, 0.0;
                     7.0, -2.0, 0.0;
                     7.0, -2.0, 0.0;

                    % Pitch Thrusters (Ty)
                   -12.0,  0.0, 0.0;
                   -12.0,  0.0, 0.0;
                     7.0,  0.0, 0.0;
                     7.0,  0.0, 0.0;
                     7.0,  0.0, 0.0;
                     7.0,  0.0, 0.0;

                    % Yaw Thrusters (Tz)
                   -12.0,  0.0, 0.3;
                   -12.0,  0.0, 0.3;
                     7.0,  0.0, 0.0;
                     7.0,  0.0, 0.0]';

    % RCS thruster direction:
    out.RCS_dir = [ % Roll Thrusters (Tx)
                    0 0  1;
                    0 0 -1;
                    0 0  1;
                    0 0 -1;

                    % Pitch Thrusters (Ty)
                    0 0  1;
                    0 0 -1;
                    0 0  1;
                    0 0  1;
                    0 0 -1;
                    0 0 -1;

                    % Yaw Thrusters (Tz)
                    0  1 0;
                    0 -1 0;
                    0  1 0;
                    0 -1 0]';
end
