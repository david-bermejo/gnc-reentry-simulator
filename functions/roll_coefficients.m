function res = roll_coefficients()
    %% Roll derivate coefficient, clean configuration
    Clb0 = [-0.00335, -0.00190, -0.00118, -0.00066, -0.00032, -0.00011, -0.00013, -0.00013, -0.00013;
            -0.00280, -0.00170, -0.00115, -0.00073, -0.00049, -0.00032, -0.00011, -0.00011, -0.00011;
            -0.00283, -0.00180, -0.00128, -0.00091, -0.00070, -0.00056, -0.00035, -0.00035, -0.00035;
            -0.00297, -0.00197, -0.00146, -0.00115, -0.00091, -0.00080, -0.00060, -0.00060, -0.00060;
            -0.00304, -0.00211, -0.00163, -0.00132, -0.00115, -0.00104, -0.00087, -0.00087, -0.00087;
            -0.00307, -0.00221, -0.00176, -0.00149, -0.00132, -0.00125, -0.00111, -0.00111, -0.00111;
            -0.00321, -0.00235, -0.00197, -0.00170, -0.00156, -0.00149, -0.00135, -0.00135, -0.00135;
            -0.00335, -0.00252, -0.00214, -0.00187, -0.00176, -0.00170, -0.00159, -0.00159, -0.00159;
            -0.00355, -0.00273, -0.00231, -0.00207, -0.00194, -0.00190, -0.00183, -0.00183, -0.00183;
            -0.00369, -0.00287, -0.00249, -0.00221, -0.00211, -0.00211, -0.00207, -0.00207, -0.00207];
    
    %% Roll coefficient increment, left elevon
    Clel = zeros(10, 9, 9);

    % Mach = 1.2
    Clel(:,:,1) = [-0.0286, -0.0243, -0.0195, -0.0111, 0.0000, 0.0056, 0.0093, 0.0112,  0.0108;
                   -0.0293, -0.0252, -0.0202, -0.0092, 0.0000, 0.0047, 0.0081, 0.0107,  0.0100;
                   -0.0293, -0.0254, -0.0177, -0.0063, 0.0000, 0.0047, 0.0081, 0.0097,  0.0090;
                   -0.0298, -0.0257, -0.0149, -0.0056, 0.0000, 0.0047, 0.0079, 0.0093,  0.0077;
                   -0.0295, -0.0234, -0.0124, -0.0056, 0.0000, 0.0047, 0.0077, 0.0084,  0.0063;
                   -0.0290, -0.0211, -0.0120, -0.0058, 0.0000, 0.0045, 0.0072, 0.0072,  0.0045;
                   -0.0285, -0.0193, -0.0122, -0.0058, 0.0000, 0.0045, 0.0065, 0.0061,  0.0024;
                   -0.0279, -0.0188, -0.0124, -0.0058, 0.0000, 0.0040, 0.0056, 0.0045, -0.0001;
                   -0.0263, -0.0193, -0.0124, -0.0058, 0.0000, 0.0038, 0.0047, 0.0024, -0.0031;
                   -0.0261, -0.0195, -0.0124, -0.0056, 0.0000, 0.0033, 0.0033, 0.0001, -0.0065];

    % Mach = 1.5
    Clel(:,:,2) =  [-0.0170, -0.0138, -0.0106, -0.0060, 0.0000, 0.0031, 0.0056, 0.0074, 0.0084;
                    -0.0170, -0.0140, -0.0113, -0.0051, 0.0000, 0.0026, 0.0054, 0.0072, 0.0079;
                    -0.0168, -0.0140, -0.0099, -0.0037, 0.0000, 0.0029, 0.0056, 0.0074, 0.0077;
                    -0.0168, -0.0142, -0.0085, -0.0033, 0.0000, 0.0031, 0.0058, 0.0074, 0.0072;
                    -0.0174, -0.0133, -0.0074, -0.0035, 0.0000, 0.0033, 0.0061, 0.0072, 0.0065;
                    -0.0181, -0.0126, -0.0074, -0.0040, 0.0000, 0.0033, 0.0058, 0.0068, 0.0056;
                    -0.0177, -0.0120, -0.0081, -0.0042, 0.0000, 0.0083, 0.0056, 0.0061, 0.0042;
                    -0.0172, -0.0122, -0.0085, -0.0042, 0.0000, 0.0033, 0.0052, 0.0049, 0.0026;
                    -0.0170, -0.0129, -0.0088, -0.0042, 0.0000, 0.0031, 0.0047, 0.0038, 0.0008;
                    -0.0174, -0.0136, -0.0090, -0.0044, 0.0000, 0.0029, 0.0038, 0.0022, 0.0012];

    % Mach = 2
    Clel(:,:,3) = [-0.0115, -0.0088, -0.0063, -0.0037, 0.0000, 0.0020, 0.0038, 0.0058, 0.0070;
                   -0.0111, -0.0085, -0.0067, -0.0031, 0.0000, 0.0020, 0.0040, 0.0061, 0.0072;
                   -0.0106, -0.0085, -0.0058, -0.0024, 0.0000, 0.0022, 0.0047, 0.0065, 0.0074;
                   -0.0104, -0.0088, -0.0053, -0.0024, 0.0000, 0.0026, 0.0052, 0.0070, 0.0077;
                   -0.0108, -0.0085, -0.0051, -0.0026, 0.0000, 0.0029, 0.0056, 0.0072, 0.0072;
                   -0.0117, -0.0085, -0.0056, -0.0031, 0.0000, 0.0031, 0.0056, 0.0070, 0.0068;
                   -0.0120, -0.0088, -0.0063, -0.0035, 0.0000, 0.0081, 0.0056, 0.0065, 0.0056;
                   -0.0124, -0.0095, -0.0069, -0.0037, 0.0000, 0.0031, 0.0052, 0.0058, 0.0045;
                   -0.0131, -0.0106, -0.0076, -0.0040, 0.0000, 0.0081, 0.0049, 0.0049, 0.0029;
                   -0.0140, -0.0113, -0.0079, -0.0040, 0.0000, 0.0029, 0.0042, 0.0038, 0.0013];

    % Mach = 3
    Clel(:,:,4) = [-0.0079, -0.0053, -0.0033, -0.0019, 0.0000, 0.0011, 0.0029, 0.0049, 0.0065;
                   -0.0072, -0.0049, -0.0037, -0.0017, 0.0000, 0.0015, 0.0036, 0.0058, 0.0072;
                   -0.0065, -0.0047, -0.0033, -0.0015, 0.0000, 0.0020, 0.0042, 0.0065, 0.0077;
                   -0.0065, -0.0051, -0.0035, -0.0019, 0.0000, 0.0024, 0.0049, 0.0070, 0.0081;
                   -0.0067, -0.0053, -0.0037, -0.0024, 0.0000, 0.0026, 0.0054, 0.0072, 0.0079;
                   -0.0076, -0.0060, -0.0044, -0.0028, 0.0000, 0.0029, 0.0058, 0.0074, 0.0077;
                   -0.0085, -0.0069, -0.0053, -0.0033, 0.0000, 0.0031, 0.0058, 0.0070, 0.0068;
                   -0.0097, -0.0081, -0.0063, -0.0035, 0.0000, 0.0033, 0.0058, 0.0065, 0.0058;
                   -0.0108, -0.0095, -0.0069, -0.0037, 0.0000, 0.0031, 0.0052, 0.0056, 0.0042;
                   -0.0122, -0.0104, -0.0074, -0.0037, 0.0000, 0.0029, 0.0049, 0.0047, 0.0026];

    % Mach = 5
    Clel(:,:,5) = [-0.0063, -0.0035, -0.0017, -0.0008, 0.0000, 0.0008, 0.0026, 0.0047, 0.0065;
                   -0.0049, -0.0028, -0.0015, -0.0008, 0.0000, 0.0013, 0.0036, 0.0058, 0.0074;
                   -0.0044, -0.0026, -0.0017, -0.0010, 0.0000, 0.0020, 0.0042, 0.0068, 0.0081;
                   -0.0047, -0.0031, -0.0024, -0.0019, 0.0000, 0.0024, 0.0052, 0.0074, 0.0086;
                   -0.0047, -0.0040, -0.0031, -0.0021, 0.0000, 0.0029, 0.0056, 0.0079, 0.0086;
                   -0.0058, -0.0051, -0.0042, -0.0028, 0.0000, 0.0031, 0.0061, 0.0079, 0.0084;
                   -0.0069, -0.0063, -0.0058, -0.0033, 0.0000, 0.0033, 0.0061, 0.0077, 0.0074;
                   -0.0085, -0.0079, -0.0063, -0.0035, 0.0000, 0.0033, 0.0058, 0.0070, 0.0065;
                   -0.0101, -0.0092, -0.0069, -0.0037, 0.0000, 0.0033, 0.0056, 0.0063, 0.0052;
                   -0.0117, -0.0104, -0.0076, -0.0040, 0.0000, 0.0031, 0.0049, 0.0049, 0.0036];

    % Mach = 10
    Clel(:,:,6) = [-0.0053, -0.0026, -0.0008, -0.0001, 0.0000, 0.0008, 0.0024, 0.0047, 0.0068;
                   -0.0042, -0.0019, -0.0005, -0.0005, 0.0000, 0.0013, 0.0036, 0.0058, 0.0079;
                   -0.0033, -0.0017, -0.0010, -0.0010, 0.0000, 0.0020, 0.0045, 0.0070, 0.0088;
                   -0.0031, -0.0021, -0.0019, -0.0015, 0.0000, 0.0024, 0.0052, 0.0077, 0.0095;
                   -0.0037, -0.0031, -0.0031, -0.0021, 0.0000, 0.0029, 0.0058, 0.0084, 0.0098;
                   -0.0047, -0.0044, -0.0042, -0.0028, 0.0000, 0.0031, 0.0063, 0.0090, 0.0101;
                   -0.0063, -0.0063, -0.0053, -0.0031, 0.0000, 0.0033, 0.0068, 0.0097, 0.0097;
                   -0.0081, -0.0079, -0.0063, -0.0035, 0.0000, 0.0036, 0.0072, 0.0095, 0.0096;
                   -0.0099, -0.0097, -0.0072, -0.0040, 0.0000, 0.0038, 0.0079, 0.0092, 0.0088;
                   -0.0122, -0.0106, -0.0079, -0.0042, 0.0000, 0.0042, 0.0090, 0.0087, 0.0058];

    % Mach = 20
    Clel(:,:,7) = [-0.0042, -0.0017, -0.0003, -0.0001, 0.0000, 0.0004, 0.0020, 0.0040, 0.0061;
                   -0.0028, -0.0010, -0.0001, -0.0001, 0.0000, 0.0011, 0.0031, 0.0054, 0.0072;
                   -0.0021, -0.0008, -0.0005, -0.0005, 0.0000, 0.0015, 0.0040, 0.0063, 0.0084;
                   -0.0019, -0.0012, -0.0012, -0.0012, 0.0000, 0.0022, 0.0049, 0.0072, 0.0088;
                   -0.0024, -0.0024, -0.0024, -0.0017, 0.0000, 0.0026, 0.0056, 0.0079, 0.0090;
                   -0.0087, -0.0037, -0.0035, -0.0024, 0.0000, 0.0029, 0.0061, 0.0081, 0.0088;
                   -0.0053, -0.0053, -0.0049, -0.0091, 0.0000, 0.0033, 0.0063, 0.0081, 0.0084;
                   -0.0069, -0.0069, -0.0058, -0.0033, 0.0000, 0.0033, 0.0063, 0.0077, 0.0074;
                   -0.0088, -0.0085, -0.0067, -0.0037, 0.0000, 0.0036, 0.0061, 0.0070, 0.0063;
                   -0.0111, -0.0099, -0.0074, -0.0040, 0.0000, 0.0033, 0.0056, 0.0061, 0.0047];

    % Mach = 30
    Clel(:,:,8) = [-0.0042, -0.0017, -0.0003, -0.0001, 0.0000, 0.0004, 0.0020, 0.0040, 0.0061;
                   -0.0028, -0.0010, -0.0001, -0.0001, 0.0000, 0.0011, 0.0031, 0.0054, 0.0072;
                   -0.0021, -0.0008, -0.0005, -0.0005, 0.0000, 0.0015, 0.0040, 0.0063, 0.0084;
                   -0.0019, -0.0012, -0.0012, -0.0012, 0.0000, 0.0022, 0.0049, 0.0072, 0.0088;
                   -0.0024, -0.0024, -0.0024, -0.0017, 0.0000, 0.0026, 0.0056, 0.0079, 0.0090;
                   -0.0087, -0.0037, -0.0035, -0.0024, 0.0000, 0.0029, 0.0061, 0.0081, 0.0088;
                   -0.0053, -0.0053, -0.0049, -0.0091, 0.0000, 0.0033, 0.0063, 0.0081, 0.0084;
                   -0.0069, -0.0069, -0.0058, -0.0033, 0.0000, 0.0033, 0.0063, 0.0077, 0.0074;
                   -0.0088, -0.0085, -0.0067, -0.0037, 0.0000, 0.0036, 0.0061, 0.0070, 0.0063;
                   -0.0111, -0.0099, -0.0074, -0.0040, 0.0000, 0.0033, 0.0056, 0.0061, 0.0047];

    % Mach = 40
    Clel(:,:,9) = [-0.0042, -0.0017, -0.0003, -0.0001, 0.0000, 0.0004, 0.0020, 0.0040, 0.0061;
                   -0.0028, -0.0010, -0.0001, -0.0001, 0.0000, 0.0011, 0.0031, 0.0054, 0.0072;
                   -0.0021, -0.0008, -0.0005, -0.0005, 0.0000, 0.0015, 0.0040, 0.0063, 0.0084;
                   -0.0019, -0.0012, -0.0012, -0.0012, 0.0000, 0.0022, 0.0049, 0.0072, 0.0088;
                   -0.0024, -0.0024, -0.0024, -0.0017, 0.0000, 0.0026, 0.0056, 0.0079, 0.0090;
                   -0.0087, -0.0037, -0.0035, -0.0024, 0.0000, 0.0029, 0.0061, 0.0081, 0.0088;
                   -0.0053, -0.0053, -0.0049, -0.0091, 0.0000, 0.0033, 0.0063, 0.0081, 0.0084;
                   -0.0069, -0.0069, -0.0058, -0.0033, 0.0000, 0.0033, 0.0063, 0.0077, 0.0074;
                   -0.0088, -0.0085, -0.0067, -0.0037, 0.0000, 0.0036, 0.0061, 0.0070, 0.0063;
                   -0.0111, -0.0099, -0.0074, -0.0040, 0.0000, 0.0033, 0.0056, 0.0061, 0.0047];

    %% Return the cell
    res = {Clb0, Clel};
end