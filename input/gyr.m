function out = gyr()
    %% Time definitions
    % Internal sample frequency:
    out.internal_freq       = 5000; % [Hz]
    % Output sample frequency:
    out.output_freq         = 200;  % [Hz]
    % Internal sample time:
    out.internal_tsample    = 1/out.internal_freq;  % [s]
    % Output sample time:
    out.output_tsample      = 1/out.output_freq;    % [s]
    
    % Output delay:
    out.output_delay = 0; % [s]

    %% Gyroscope model error properties
    % Maximum angular velocity:
    out.max_ang_vel = 250; % [deg/s]
    
    % Scale factor:
    out.sf          = 10 * 1e-6 * 0;    % [ppm]
    % Asymmetric scale factor:
    out.sf_asym     = 0 * 1e-6;         % [ppm]
    % Nonlinear scale factor:
    out.sf_nl       = 20 * 1e-6 * 0;    % [ppm]

    % Groscope bias:
    bias_rep        = 0.01; % [deg/h]
    out.bias        = randn(3,1) * (bias_rep/3600); % [deg/s]

    % Angular Random Walk:
    N_ARW           = 0.005;            % [deg/sqrt(h)]
    out.arw_mu      = 0;
    out.arw_var     = (N_ARW/60)^2;     % [(deg/sqrt(s))^2]
    out.arw_seed    = randi(100000,3,1);

    % Bias Instability:
    B_ins           = 0.005;            % [deg/h]
    out.ins_Tc      = 300;              % [s]
    out.ins_var     = (B_ins/3600)^2;   % [(deg/s)^2]
    out.ins_seed    = randi(100000,3,1);

    % Rate Random Walk:
    N_RRW           = 0.005;            % [deg/h^(3/2)]
    out.rrw_mu      = 0;
    out.rrw_var     = (N_RRW/60^3)^2;   % [(deg/s^2)^2]
    out.rrw_seed    = randi(100000,3,1);

    % Quatization resolution:
    out.quant_res   = out.max_ang_vel/2^32; % [deg/bit]

    % Mounting error DCM:
    out.dcm_mnt_err = [1 0 0; 0 1 0; 0 0 1];
    % Gyr to Body frame DCM
    out.dcm_gyr2bod = [1 0 0; 0 1 0; 0 0 1];
    % Non-orthogonality error DCM:
    out.non_orth    = [1 0 0; 0 1 0; 0 0 1];
end