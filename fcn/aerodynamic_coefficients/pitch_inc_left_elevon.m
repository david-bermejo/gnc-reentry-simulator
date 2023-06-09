function Cmel = pitch_inc_left_elevon(AoA, delta_el, M)
    AoA_base = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]; % deg
    delta_el_base = [-40, -30, -20, -10, 0, 10, 20, 30, 40]; % deg
    M_base = [1.2, 1.5, 2, 3, 5, 10, 20];
    
    % Mach = 1.2
    Cmel_base = [0.0263, 0.0234, 0.0193, 0.0114, 0.0000, -0.0067, -0.0120, -0.0166, -0.0193;
                 0.0278, 0.0246, 0.0203, 0.0094, 0.0000, -0.0058, -0.0113, -0.0156, -0.0181;
                 0.0283, 0.0251, 0.0179, 0.0068, 0.0000, -0.0062, -0.0115, -0.0157, -0.0178;
                 0.0290, 0.0258, 0.0155, 0.0061, 0.0000, -0.0063, -0.0116, -0.0156, -0.0171;
                 0.0304, 0.0239, 0.0135, 0.0063, 0.0000, -0.0065, -0.0118, -0.0152, -0.0159;
                 0.0314, 0.0224, 0.0133, 0.0067, 0.0000, -0.0065, -0.0116, -0.0145, -0.0145;
                 0.0304, 0.0208, 0.0137, 0.0068, 0.0000, -0.0063, -0.0113, -0.0135, -0.0130;
                 0.0294, 0.0208, 0.0142, 0.0070, 0.0000, -0.0062, -0.0104, -0.0123, -0.0103;
                 0.0285, 0.0215, 0.0145, 0.0072, 0.0000, -0.0062, -0.0097, -0.0104, -0.0074;
                 0.0287, 0.0220, 0.0147, 0.0068, 0.0000, -0.0056, -0.0086, -0.0084, -0.0039];

    % Mach = 1.5
    Cmel_base(:,:,2) = [0.0155, 0.0131, 0.0104, 0.0060, 0.0000, -0.0039, -0.0072, -0.0106, -0.0132;
                        0.0159, 0.0135, 0.0111, 0.0051, 0.0000, -0.0036, -0.0072, -0.0106, -0.0132;
                        0.0159, 0.0137, 0.0097, 0.0036, 0.0000, -0.0039, -0.0079, -0.0111, -0.0133;
                        0.0162, 0.0142, 0.0087, 0.0034, 0.0000, -0.0043, -0.0082, -0.0116, -0.0135;
                        0.0169, 0.0135, 0.0080, 0.0039, 0.0000, -0.0045, -0.0086, -0.0116, -0.0132;
                        0.0181, 0.0131, 0.0082, 0.0043, 0.0000, -0.0046, -0.0087, -0.0115, -0.0123;
                        0.0176, 0.0128, 0.0087, 0.0046, 0.0000, -0.0048, -0.0087, -0.0111, -0.0113;
                        0.0176, 0.0131, 0.0094, 0.0048, 0.0000, -0.0048, -0.0084, -0.0103, -0.0099;
                        0.0181, 0.0142, 0.0099, 0.0051, 0.0000, -0.0046, -0.0080, -0.0092, -0.0080;
                        0.0188, 0.0150, 0.0104, 0.0049, 0.0000, -0.0045, -0.0074, -0.0077, -0.0060];

    % Mach = 2
    Cmel_base(:,:,3) = [0.0104, 0.0082, 0.0061, 0.0034, 0.0000, -0.0022, -0.0051, -0.0080, -0.0106;
                        0.0102, 0.0080, 0.0063, 0.0029, 0.0000, -0.0024, -0.0056, -0.0086, -0.0111;
                        0.0099, 0.0082, 0.0060, 0.0024, 0.0000, -0.0029, -0.0062, -0.0092, -0.0118;
                        0.0101, 0.0087, 0.0056, 0.0024, 0.0000, -0.0033, -0.0070, -0.0097, -0.0121;
                        0.0108, 0.0085, 0.0055, 0.0029, 0.0000, -0.0036, -0.0074, -0.0104, -0.0123;
                        0.0116, 0.0089, 0.0060, 0.0032, 0.0000, -0.0038, -0.0077, -0.0106, -0.0120;
                        0.0121, 0.0090, 0.0067, 0.0037, 0.0000, -0.0043, -0.0079, -0.0103, -0.0115;
                        0.0126, 0.0101, 0.0075, 0.0041, 0.0000, -0.0043, -0.0079, -0.0099, -0.0103;
                        0.0135, 0.0111, 0.0084, 0.0041, 0.0000, -0.0043, -0.0075, -0.0091, -0.0089;
                        0.0147, 0.0121, 0.0085, 0.0046, 0.0000, -0.0041, -0.0070, -0.0084, -0.0070];

    % Mach = 3
    Cmel_base(:,:,4) = [0.0072, 0.0048, 0.0031, 0.0015, 0.0000, -0.0015, -0.0039, -0.0067, -0.0092;
                        0.0065, 0.0044, 0.0031, 0.0014, 0.0000, -0.0019, -0.0046, -0.0077, -0.0103;
                        0.0060, 0.0044, 0.0031, 0.0014, 0.0000, -0.0024, -0.0055, -0.0086, -0.0113;
                        0.0060, 0.0049, 0.0032, 0.0017, 0.0000, -0.0029, -0.0063, -0.0097, -0.0120;
                        0.0065, 0.0053, 0.0037, 0.0022, 0.0000, -0.0036, -0.0072, -0.0103, -0.0123;
                        0.0075, 0.0063, 0.0046, 0.0027, 0.0000, -0.0038, -0.0075, -0.0106, -0.0123;
                        0.0085, 0.0070, 0.0055, 0.0032, 0.0000, -0.0041, -0.0079, -0.0104, -0.0118;
                        0.0097, 0.0082, 0.0065, 0.0036, 0.0000, -0.0043, -0.0079, -0.0103, -0.0109;
                        0.0111, 0.0097, 0.0075, 0.0039, 0.0000, -0.0043, -0.0075, -0.0094, -0.0094;
                        0.0125, 0.0109, 0.0078, 0.0041, 0.0000, -0.0043, -0.0072, -0.0087, -0.0084];

    % Mach = 5
    Cmel_base(:,:,5) = [0.0055, 0.0031, 0.0015, 0.0007, 0.0000, -0.0012, -0.0033, -0.0062, -0.0089;
                        0.0044, 0.0026, 0.0014, 0.0007, 0.0000, -0.0017, -0.0043, -0.0077, -0.0104;
                        0.0039, 0.0024, 0.0015, 0.0008, 0.0000, -0.0022, -0.0055, -0.0087, -0.0116;
                        0.0041, 0.0031, 0.0024, 0.0015, 0.0000, -0.0029, -0.0063, -0.0097, -0.0123;
                        0.0046, 0.0039, 0.0032, 0.0022, 0.0000, -0.0033, -0.0070, -0.0104, -0.0127;
                        0.0056, 0.0049, 0.0043, 0.0031, 0.0000, -0.0038, -0.0075, -0.0108, -0.0127;
                        0.0073, 0.0067, 0.0055, 0.0034, 0.0000, -0.0041, -0.0077, -0.0106, -0.0121;
                        0.0087, 0.0080, 0.0065, 0.0037, 0.0000, -0.0041, -0.0079, -0.0104, -0.0113;
                        0.0106, 0.0096, 0.0075, 0.0043, 0.0000, -0.0043, -0.0077, -0.0097, -0.0099;
                        0.0123, 0.0111, 0.0082, 0.0044, 0.0000, -0.0041, -0.0072, -0.0087, -0.0086];

    % Mach = 10
    Cmel_base(:,:,6) = [0.0044, 0.0019, 0.0003, 0.0000, 0.0000, -0.0012, -0.0036, -0.0063, -0.0091;
                        0.0034, 0.0014, 0.0003, 0.0002, 0.0000, -0.0019, -0.0045, -0.0077, -0.0108;
                        0.0029, 0.0014, 0.0008, 0.0008, 0.0000, -0.0022, -0.0055, -0.0087, -0.0121;
                        0.0029, 0.0019, 0.0017, 0.0015, 0.0000, -0.0029, -0.0063, -0.0101, -0.0135;
                        0.0034, 0.0029, 0.0029, 0.0019, 0.0000, -0.0034, -0.0070, -0.0111, -0.0149;
                        0.0046, 0.0044, 0.0043, 0.0026, 0.0000, -0.0038, -0.0079, -0.0120, -0.0162;
                        0.0061, 0.0061, 0.0051, 0.0031, 0.0000, -0.0043, -0.0086, -0.0133, -0.0176;
                        0.0080, 0.0077, 0.0061, 0.0034, 0.0000, -0.0045, -0.0092, -0.0149, -0.0174;
                        0.0101, 0.0094, 0.0073, 0.0039, 0.0000, -0.0050, -0.0104, -0.0159, -0.0152;
                        0.0121, 0.0108, 0.0082, 0.0044, 0.0000, -0.0055, -0.0121, -0.0157, -0.0118];

    % Mach = 20
    Cmel_base(:,:,7) = [0.0036, 0.0014, 0.0000, 0.0000, 0.0000, -0.0005, -0.0024, -0.0053, -0.0080;
                        0.0024, 0.0005, 0.0000, 0.0000, 0.0000, -0.0015, -0.0038, -0.0070, -0.0101;
                        0.0017, 0.0003, 0.0002, 0.0002, 0.0000, -0.0022, -0.0050, -0.0084, -0.0113;
                        0.0017, 0.0012, 0.0012, 0.0012, 0.0000, -0.0026, -0.0060, -0.0094, -0.0123;
                        0.0022, 0.0022, 0.0022, 0.0017, 0.0000, -0.0033, -0.0068, -0.0104, -0.0130;
                        0.0036, 0.0036, 0.0036, 0.0022, 0.0000, -0.0036, -0.0074, -0.0109, -0.0132;
                        0.0053, 0.0053, 0.0049, 0.0027, 0.0000, -0.0041, -0.0079, -0.0111, -0.0128;
                        0.0070, 0.0070, 0.0058, 0.0034, 0.0000, -0.0043, -0.0082, -0.0109, -0.0121;
                        0.0087, 0.0085, 0.0068, 0.0037, 0.0000, -0.0045, -0.0080, -0.0104, -0.0111;
                        0.0109, 0.0099, 0.0077, 0.0041, 0.0000, -0.0045, -0.0077, -0.0097, -0.0097];

    Cmel = interp3(delta_el_base, AoA_base, M_base, Cmel_base, delta_el, AoA, M, 'spline');
end