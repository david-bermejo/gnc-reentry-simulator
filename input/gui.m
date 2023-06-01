function out = gui()
    %% Time definitions
    % Execution frequency:
    out.freq = 20;            % [Hz]
    % Execution time interval:
    out.tsample = 1/out.freq; % [s]
    % Output delay:
    out.delay = 0.005; % [s]

    %% Guidance Linear Quadratic Regulator configuration
    % Maximum radial position deviation:
    R_max      = 100;             % [m]
    % Maximum latitude deviation:
    delta_max  = deg2rad(1.0);    % [rad]
    % Maximum velocity deviation (magnitude):
    V_max      = 5;               % [m/s]
    % Maximum glideslope deviation:
    gamma_max  = deg2rad(0.1);   % [rad]
    % Maximum heading deviation:
    chi_max    = deg2rad(1.5);    % [rad]
    % Maximum commanded angle of attack deviation:
    AoA_max  = deg2rad(0.5);        % [rad]
    % Maximum commanded bank angle deviation:
    sigma_max  = deg2rad(0.5); % [rad]

    % State cost matrix:
    out.Q      = zeros(3);
    out.Q(1,1) = 1/R_max^2;
    %out.Q(2,2) = 1/delta_max^2;
    out.Q(2,2) = 1/V_max^2;
    out.Q(3,3) = 1/gamma_max^2;
    %out.Q(5,5) = 1/chi_max^2;

    % Control cost matrix:
    out.R      = zeros(2);
    out.R(1,1) = 1/AoA_max^2;
    out.R(2,2) = 1/sigma_max^2;

    %% Reference Trajectory (offline)
    load('trajectory-v33.mat', 'ts', 'y', 'u', 'ef');
    % Reference time:
    out.t_gui = ts;
    % Reference state:
    out.y_gui = y;
    % Reference control input (bank & body flap):
    out.u_gui = u;

    % Target Longitude:
    out.tauf        = deg2rad(77.45); % [rad]
    % Target Latitude:
    out.deltaf      = deg2rad(18.44); % [rad]
    % Maximum Heading deviation interpolant:
    out.d_chi0    = deg2rad(4);       % [rad]
    out.d_chif    = deg2rad(2.5);       % [rad]
    out.M0        = 40;               % [-]
    out.Mf        = 5;                % [-]
    % Maximum bank rate:
    out.dbank_max   = deg2rad(15);    % [rad/s]
    % TAEM energy transition:
    out.ef          = ef;             % [J/kg]

    %% Workspace cleanup
    vars = {'ts', 'y', 'u'};
    clear('base', vars{:});
end