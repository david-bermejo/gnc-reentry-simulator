function out = gui()
    R_max       = 100;              % [m]
    tau_max     = deg2rad(0.5);     % [rad]
    delta_max   = deg2rad(0.5);     % [rad]
    V_max       = 50;               % [m/s]
    gamma_max   = deg2rad(0.5);    % [rad]
    chi_max     = deg2rad(2.5);     % [rad]
    sigma_max   = deg2rad(2.5);       % [rad]

    Q = zeros(3);
    Q(1,1) = 1/R_max^2;
    %Q(2,2) = 1/tau_max^2;
    %Q(2,2) = 1/delta_max^2;
    Q(2,2) = 1/V_max^2;
    Q(3,3) = 1/gamma_max^2;
    %Q(5,5) = 1/chi_max^2;

    out.Q = Q;
    out.R = [1/sigma_max^2];

    % Guidance Algorithm frequency [Hz]
    out.freq = 0.5;
    out.tsample = 1/out.freq;

    % Delay [s]
    out.delay = 0.05;

    % Reference Guidance
    load('trajectory.mat', 'ts', 'y', 'energy', 'u');
    %elems = energy < 1.0;
    %out.energy_gui = flip(min(energy(elems), 1.0));
    %out.u_gui = flip(u(:,elems), 2);
    out.t_gui = ts;
    out.y_gui = y;
    out.u_gui = u;

    vars = {'ts', 'y', 'energy', 'u'};
    clear('base', vars{:});
end