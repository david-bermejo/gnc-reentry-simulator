function out = rad(noiseFlag)
    %% Constants
    % Planet radius:
    Rp = 3396.2 * 1e3; % [m]

    %% Time characteristics
    % Sample frequency:
    out.freq = 10;          % [Hz]
    % Sample time:
    out.tsamp = 1/out.freq; % [s]

    % Beacons longitude and latitude (Mars surface):
    bc_ang_pos_ecef = deg2rad([74, 16; ...
                               79, 17; ...
                               80, 20; ...
                               75, 21])'; % [deg]

    % Number of beacons considered:
    out.N_beacons = size(bc_ang_pos_ecef, 2);

    % Beacon positions (ECEF frame):
    out.bc_pos_ecef = [Rp.*cos(bc_ang_pos_ecef(1,:)).*cos(bc_ang_pos_ecef(2,:));
                       Rp.*sin(bc_ang_pos_ecef(1,:)).*cos(bc_ang_pos_ecef(2,:));
                       Rp.*sin(bc_ang_pos_ecef(2,:))]; % [m]

    %% Noise definitions
    % Pseudorange variance:
    out.rand_rang_var = 10^2;
    % Pseudorange seed (RNG):
    out.rand_rang_seed = randi(10000,1);

    % Pseudorange rate variance:
    out.rand_dopl_var = 0.1^2;
    % Pseudorange rate seed (RNG):
    out.rand_dopl_seed = randi(10000,1);
end