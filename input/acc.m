function out = acc(noiseFlag)
    %% ACC time constants
    % Accelerometer sampling frequency [Hz]
    out.internal_freq       = 5000;
    out.output_freq         = 200;
    out.internal_tsample    = 1/out.internal_freq;
    out.output_tsample      = 1/out.output_freq;
    
    % Output delay [s]
    %out.output_delay = out.output_tsample;
    out.output_delay = 0;

    %% ACC model error properties
    % Gravity Reference
    out.g0 = 9.80665; % [m/s^2]

    % Max acceleration [g's]
    out.max_acc = 10 * out.g0;

    % Scale factors [ppm]
    out.sf          = 30 * 1e-6 * 0;
    out.sf_asym     = 0 * 1e-6;
    out.sf_nl       = 150 * 1e-6 * 0;

    % Bias [m/s^2]
    bias_rep        = 0.1*1e-3*out.g0; % [ug]
    out.bias        = randn(3,1) * bias_rep;

    % Velocity Random Walk [m/s^(3/2)]
    N_VRW           = 0.1; % [ug/sqrt(Hz)]
    out.vrw_mu      = 0;
    out.vrw_var     = (N_VRW*1e-3*out.g0)^2;
    out.vrw_seed    = randi(100000,3,1);

    % Bias Instability [m/s^2]
    B_ins           = 0.01; % [ug]
    out.ins_Tc      = 300; % [s]
    out.ins_var     = (B_ins*1e-3*out.g0)^2;
    out.ins_seed    = randi(100000,3,1);

    % Acceleration Random Walk [m/s^(5/2)]
    N_ARW           = 0.05; % [ug/sqrt(hr)]
    out.arw_mu      = 0;
    out.arw_var     = (N_ARW*1e-3*out.g0/60)^2;
    out.arw_seed    = randi(100000,3,1);

    % Quatization resolution
    out.quant_res   = out.max_acc/2^32;

    out.dcm_mnt_err = [1 0 0; 0 1 0; 0 0 1];
    out.dcm_acc2bod = [1 0 0; 0 1 0; 0 0 1];
    out.non_orth    = [1 0 0; 0 1 0; 0 0 1];
end