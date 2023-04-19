function Csel = side_inc_left_elevon(AoA, delta_el, M)
    AoA_base = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]; % deg
    delta_el_base = [-40, -30, -20, -10, 0, 10, 20, 30, 40]; % deg
    M_base = [1.2, 1.5, 2, 3, 5, 10, 20];

    % Mach = 1.2
    Csel_base = [-0.0053, -0.0074, -0.0081, -0.0057, 0.0000, 0.0046, 0.0098, 0.0156, 0.0214;
                 -0.0070, -0.0091, -0.0091, -0.0053, 0.0000, 0.0046, 0.0101, 0.0160, 0.0221;
                 -0.0084, -0.0101, -0.0087, -0.0046, 0.0000, 0.0050, 0.0108, 0.0173, 0.0228;
                 -0.0101, -0.0111, -0.0087, -0.0043, 0.0000, 0.0053, 0.0118, 0.0186, 0.0238;
                 -0.0118, -0.0118, -0.0087, -0.0050, 0.0000, 0.0057, 0.0125, 0.0190, 0.0245;
                 -0.0132, -0.0125, -0.0094, -0.0053, 0.0000, 0.0063, 0.0129, 0.0194, 0.0245;
                 -0.0149, -0.0132, -0.0105, -0.0057, 0.0000, 0.0063, 0.0135, 0.0197, 0.0242;
                 -0.0166, -0.0146, -0.0111, -0.0063, 0.0000, 0.0067, 0.0139, 0.0197, 0.0231;
                 -0.0183, -0.0159, -0.0118, -0.0067, 0.0000, 0.0070, 0.0139, 0.0190, 0.0221;
                 -0.0201, -0.0173, -0.0129, -0.0070, 0.0000, 0.0070, 0.0135, 0.0183, 0.0201];

    % Mach = 1.5
    Csel_base(:,:,2) = [-0.0029, -0.0046, -0.0046, -0.0033, 0.0000, 0.0026, 0.0060, 0.0098, 0.0146;
                        -0.0039, -0.0050, -0.0050, -0.0029, 0.0000, 0.0026, 0.0063, 0.0108, 0.0153;
                        -0.0046, -0.0057, -0.0050, -0.0026, 0.0000, 0.0029, 0.0070, 0.0118, 0.0163;
                        -0.0057, -0.0063, -0.0053, -0.0026, 0.0000, 0.0036, 0.0077, 0.0125, 0.0173;
                        -0.0067, -0.0067, -0.0053, -0.0029, 0.0000, 0.0039, 0.0087, 0.0135, 0.0180;
                        -0.0081, -0.0077, -0.0057, -0.0036, 0.0000, 0.0043, 0.0091, 0.0139, 0.0183;
                        -0.0091, -0.0081, -0.0063, -0.0039, 0.0000, 0.0046, 0.0098, 0.0142, 0.0183;
                        -0.0105, -0.0094, -0.0074, -0.0043, 0.0000, 0.0050, 0.0098, 0.0146, 0.0180;
                        -0.0118, -0.0105, -0.0081, -0.0046, 0.0000, 0.0050, 0.0101, 0.0142, 0.0170;
                        -0.0129, -0.0115, -0.0087, -0.0050, 0.0000, 0.0050, 0.0101, 0.0139, 0.0159];

    % Mach = 2
    Csel_base(:,:,3) = [-0.0019, -0.0026, -0.0026, -0.0019, 0.0000, 0.0015, 0.0039, 0.0070, 0.0111;
                        -0.0022, -0.0029, -0.0029, -0.0019, 0.0000, 0.0015, 0.0046, 0.0081, 0.0122;
                        -0.0029, -0.0033, -0.0033, -0.0015, 0.0000, 0.0022, 0.0053, 0.0094, 0.0142;
                        -0.0036, -0.0039, -0.0033, -0.0015, 0.0000, 0.0026, 0.0063, 0.0105, 0.0146;
                        -0.0046, -0.0046, -0.0036, -0.0022, 0.0000, 0.0029, 0.0070, 0.0111, 0.0153;
                        -0.0053, -0.0053, -0.0039, -0.0026, 0.0000, 0.0036, 0.0077, 0.0118, 0.0156;
                        -0.0063, -0.0060, -0.0050, -0.0029, 0.0000, 0.0036, 0.0081, 0.0122, 0.0156;
                        -0.0074, -0.0067, -0.0057, -0.0036, 0.0000, 0.0039, 0.0084, 0.0125, 0.0156;
                        -0.0087, -0.0077, -0.0063, -0.0036, 0.0000, 0.0043, 0.0084, 0.0122, 0.0153;
                        -0.0098, -0.0087, -0.0070, -0.0043, 0.0000, 0.0043, 0.0084, 0.0122, 0.0142];

    % Mach = 3
    Csel_base(:,:,4) = [-0.0012, -0.0012, -0.0012, -0.0012, 0.0000, 0.0009, 0.0029, 0.0060, 0.0098;
                        -0.0015, -0.0015, -0.0015, -0.0012, 0.0000, 0.0012, 0.0036, 0.0070, 0.0111;
                        -0.0019, -0.0019, -0.0019, -0.0009, 0.0000, 0.0019, 0.0046, 0.0084, 0.0122;
                        -0.0022, -0.0022, -0.0022, -0.0012, 0.0000, 0.0022, 0.0053, 0.0094, 0.0132;
                        -0.0029, -0.0029, -0.0026, -0.0015, 0.0000, 0.0026, 0.0063, 0.0101, 0.0142;
                        -0.0039, -0.0039, -0.0033, -0.0022, 0.0000, 0.0029, 0.0067, 0.0111, 0.0146;
                        -0.0046, -0.0046, -0.0039, -0.0026, 0.0000, 0.0033, 0.0074, 0.0115, 0.0149;
                        -0.0057, -0.0057, -0.0046, -0.0029, 0.0000, 0.0036, 0.0077, 0.0118, 0.0146;
                        -0.0070, -0.0067, -0.0053, -0.0033, 0.0000, 0.0039, 0.0081, 0.0118, 0.0142;
                        -0.0081, -0.0074, -0.0060, -0.0036, 0.0000, 0.0039, 0.0084, 0.0115, 0.0139];

    % Mach = 5
    Csel_base(:,:,5) = [-0.0005, -0.0005, -0.0005, -0.0005, 0.0000, 0.0005, 0.0026, 0.0053, 0.0091;
                        -0.0005, -0.0005, -0.0005, -0.0005, 0.0000, 0.0009, 0.0033, 0.0063, 0.0105;
                        -0.0009, -0.0009, -0.0009, -0.0009, 0.0000, 0.0015, 0.0043, 0.0077, 0.0118;
                        -0.0015, -0.0015, -0.0015, -0.0012, 0.0000, 0.0019, 0.0050, 0.0087, 0.0129;
                        -0.0022, -0.0022, -0.0022, -0.0015, 0.0000, 0.0022, 0.0060, 0.0098, 0.0139;
                        -0.0029, -0.0029, -0.0029, -0.0019, 0.0000, 0.0029, 0.0067, 0.0105, 0.0142;
                        -0.0039, -0.0039, -0.0036, -0.0022, 0.0000, 0.0033, 0.0070, 0.0111, 0.0146;
                        -0.0050, -0.0050, -0.0043, -0.0026, 0.0000, 0.0036, 0.0074, 0.0115, 0.0142;
                        -0.0060, -0.0060, -0.0050, -0.0033, 0.0000, 0.0036, 0.0077, 0.0115, 0.0142;
                        -0.0074, -0.0070, -0.0060, -0.0036, 0.0000, 0.0039, 0.0077, 0.0111, 0.0135];

    % Mach = 10
    Csel_base(:,:,6) = [ 0.0000,  0.0000,  0.0000,  0.0000, 0.0000, 0.0009, 0.0022, 0.0050, 0.0087;
                         0.0000,  0.0000,  0.0000,  0.0000, 0.0000, 0.0012, 0.0033, 0.0063, 0.0105;
                        -0.0005, -0.0005, -0.0005, -0.0005, 0.0000, 0.0015, 0.0043, 0.0077, 0.0122;
                        -0.0012, -0.0012, -0.0012, -0.0012, 0.0000, 0.0019, 0.0050, 0.0087, 0.0139;
                        -0.0019, -0.0019, -0.0019, -0.0015, 0.0000, 0.0022, 0.0060, 0.0101, 0.0156;
                        -0.0029, -0.0029, -0.0029, -0.0019, 0.0000, 0.0029, 0.0067, 0.0115, 0.0177;
                        -0.0036, -0.0036, -0.0034, -0.0022, 0.0000, 0.0033, 0.0074, 0.0129, 0.0194;
                        -0.0046, -0.0046, -0.0041, -0.0026, 0.0000, 0.0036, 0.0084, 0.0146, 0.0197;
                        -0.0060, -0.0060, -0.0048, -0.0033, 0.0000, 0.0039, 0.0094, 0.0159, 0.0187;
                        -0.0074, -0.0070, -0.0055, -0.0036, 0.0000, 0.0046, 0.0111, 0.0166, 0.0159];

    % Mach = 20
    Csel_base(:,:,7) = [ 0.0000,  0.0000,  0.0000,  0.0000, 0.0000, 0.0002, 0.0015, 0.0043, 0.0077;
                         0.0000,  0.0000,  0.0000,  0.0000, 0.0000, 0.0009, 0.0026, 0.0057, 0.0094;
                        -0.0005, -0.0005, -0.0005, -0.0005, 0.0000, 0.0012, 0.0036, 0.0070, 0.0111;
                        -0.0009, -0.0009, -0.0009, -0.0009, 0.0000, 0.0015, 0.0046, 0.0081, 0.0122;
                        -0.0015, -0.0015, -0.0015, -0.0012, 0.0000, 0.0022, 0.0057, 0.0091, 0.0135;
                        -0.0022, -0.0022, -0.0022, -0.0015, 0.0000, 0.0026, 0.0063, 0.0101, 0.0142;
                        -0.0033, -0.0033, -0.0029, -0.0019, 0.0000, 0.0029, 0.0067, 0.0108, 0.0146;
                        -0.0043, -0.0043, -0.0038, -0.0026, 0.0000, 0.0033, 0.0074, 0.0111, 0.0146;
                        -0.0053, -0.0053, -0.0043, -0.0029, 0.0000, 0.0036, 0.0077, 0.0115, 0.0142;
                        -0.0067, -0.0063, -0.0053, -0.0033, 0.0000, 0.0036, 0.0077, 0.0115, 0.0135];

    Csel = interp3(delta_el_base, AoA_base, M_base, Csel_base, delta_el, AoA, M, 'spline');
end