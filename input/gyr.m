function out = gyr()
    %% Gyro time constants
    % Gyroscope sampling frequency [Hz]
    out.internal_freq       = 7000;
    out.output_freq         = 20;
    out.internal_tsample    = 1/out.internal_freq;
    out.output_tsample      = 1/out.output_freq;
    
    % Output delay [s]
    out.output_delay = out.output_tsample;

    %% Gyro model error properties
    % Scale factors [ppm]
    out.sf          = 10 * 1e-6 * 0;
    out.sf_asym     = 10 * 1e-6 * 0;
    out.sf_nl       = 10 * 1e-6 * 0;

    % Bias [rad/s]
    out.bias        = 0;

    out.quant_res   = 1/2^32;

    % Max angular velocity [rad/s]
    out.max_ang_vel = deg2rad(250);

    out.dcm_mnt_err = [1 0 0; 0 1 0; 0 0 1];
    out.dcm_gyr2bod = [1 0 0; 0 1 0; 0 0 1];
    out.non_orth    = [1 0 0; 0 1 0; 0 0 1];

    out.white_noise_mu      = 0;
    out.white_noise_var     = 0;
    out.white_noise_seed    = 0;

    out.bias_inst_mu        = 0;
    out.bias_inst_var       = 0;
    out.bias_inst_seed      = 0;

    out.ARW_mu              = 0;
    out.ARW_var             = 0;
    out.ARW_seed            = 0;
end