function Cnel = yaw_inc_left_elevon(AoA, delta_el, M)
    AoA_base = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]; % deg
    delta_el_base = [-40, -30, -20, -10, 0, 10, 20, 30, 40]; % deg
    M_base = [1.2, 1.5, 2, 3, 5, 10, 20];

    % Mach = 1.2
    Cnel_base = [-0.0110, -0.0027, 0.0024, 0.0041, 0.0000, -0.0068, -0.0154, -0.0260, -0.0377;
                 -0.0072,  0.0003, 0.0048, 0.0048, 0.0000, -0.0072, -0.0164, -0.0277, -0.0397;
                 -0.0038,  0.0031, 0.0065, 0.0055, 0.0000, -0.0075, -0.0178, -0.0295, -0.0418;
                  0.0000,  0.0058, 0.0086, 0.0062, 0.0000, -0.0086, -0.0192, -0.0315, -0.0435;
                  0.0031,  0.0089, 0.0103, 0.0068, 0.0000, -0.0092, -0.0205, -0.0329, -0.0452;
                  0.0068,  0.0123, 0.0120, 0.0075, 0.0000, -0.0099, -0.0219, -0.0346, -0.0466;
                  0.0116,  0.0158, 0.0137, 0.0082, 0.0000, -0.0110, -0.0233, -0.0360, -0.0473;
                  0.0161,  0.0182, 0.0151, 0.0089, 0.0000, -0.0113, -0.0240, -0.0370, -0.0473;
                  0.0205,  0.0202, 0.0168, 0.0096, 0.0000, -0.0120, -0.0250, -0.0373, -0.0462;
                  0.0243,  0.0226, 0.0182, 0.0106, 0.0000, -0.0123, -0.0257, -0.0373, -0.0445];

    % Mach = 1.5
    Cnel_base(:,:,2) = [-0.0062, -0.0014, 0.0014, 0.0024, 0.0000, -0.0038, -0.0086, -0.0154, -0.0233;
                        -0.0041,  0.0000, 0.0024, 0.0024, 0.0000, -0.0041, -0.0096, -0.0168, -0.0250;
                        -0.0021,  0.0017, 0.0038, 0.0031, 0.0000, -0.0045, -0.0106, -0.0185, -0.0271;
                         0.0000,  0.0031, 0.0048, 0.0031, 0.0000, -0.0051, -0.0116, -0.0202, -0.0288;
                         0.0017,  0.0048, 0.0055, 0.0038, 0.0000, -0.0055, -0.0130, -0.0216, -0.0301;
                         0.0038,  0.0068, 0.0065, 0.0041, 0.0000, -0.0062, -0.0140, -0.0229, -0.0305;
                         0.0062,  0.0086, 0.0075, 0.0048, 0.0000, -0.0068, -0.0151, -0.0240, -0.0318;
                         0.0089,  0.0103, 0.0086, 0.0051, 0.0000, -0.0075, -0.0161, -0.0243, -0.0325;
                         0.0116,  0.0116, 0.0099, 0.0058, 0.0000, -0.0079, -0.0168, -0.0253, -0.0325;
                         0.0137,  0.0130, 0.0110, 0.0065, 0.0000, -0.0082, -0.0171, -0.0257, -0.0318];

    % Mach = 2
    Cnel_base(:,:,3) = [-0.0041, -0.0010, 0.0010, 0.0014, 0.0000, -0.0021, -0.0055, -0.0106, -0.0164;
                        -0.0024,  0.0000, 0.0014, 0.0014, 0.0000, -0.0024, -0.0062, -0.0120, -0.0185;
                        -0.0014,  0.0007, 0.0021, 0.0017, 0.0000, -0.0027, -0.0075, -0.0134, -0.0202;
                         0.0000,  0.0021, 0.0027, 0.0021, 0.0000, -0.0034, -0.0086, -0.0147, -0.0219;
                         0.0010,  0.0027, 0.0031, 0.0024, 0.0000, -0.0038, -0.0096, -0.0161, -0.0236;
                         0.0024,  0.0041, 0.0041, 0.0027, 0.0000, -0.0045, -0.0106, -0.0175, -0.0247;
                         0.0038,  0.0055, 0.0048, 0.0031, 0.0000, -0.0051, -0.0116, -0.0188, -0.0257;
                         0.0058,  0.0065, 0.0058, 0.0038, 0.0000, -0.0055, -0.0123, -0.0199, -0.0264;
                         0.0075,  0.0075, 0.0065, 0.0041, 0.0000, -0.0062, -0.0130, -0.0202, -0.0267;
                         0.0086,  0.0086, 0.0075, 0.0048, 0.0000, -0.0065, -0.0137, -0.0209, -0.0264];

    % Mach = 3
    Cnel_base(:,:,4) = [-0.0027, -0.0003, 0.0007, 0.0007, 0.0000, -0.0010, -0.0034, -0.0075, -0.0130;
                        -0.0021,  0.0000, 0.0007, 0.0007, 0.0000, -0.0017, -0.0045, -0.0092, -0.0154;
                        -0.0007,  0.0007, 0.0010, 0.0010, 0.0000, -0.0021, -0.0055, -0.0103, -0.0171;
                         0.0000,  0.0010, 0.0014, 0.0014, 0.0000, -0.0024, -0.0065, -0.0127, -0.0192;
                         0.0007,  0.0017, 0.0017, 0.0017, 0.0000, -0.0031, -0.0079, -0.0140, -0.0209;
                         0.0017,  0.0024, 0.0024, 0.0017, 0.0000, -0.0034, -0.0089, -0.0154, -0.0219;
                         0.0027,  0.0034, 0.0034, 0.0024, 0.0000, -0.0041, -0.0099, -0.0164, -0.0229;
                         0.0038,  0.0041, 0.0041, 0.0027, 0.0000, -0.0048, -0.0106, -0.0175, -0.0236;
                         0.0051,  0.0051, 0.0048, 0.0034, 0.0000, -0.0051, -0.0116, -0.0182, -0.0240;
                         0.0062,  0.0062, 0.0058, 0.0038, 0.0000, -0.0055, -0.0123, -0.0188, -0.0240];

    % Mach = 5
    Cnel_base(:,:,5) = [-0.0024, -0.0003, 0.0000, 0.0000, 0.0000, -0.0007, -0.0027, -0.0065, -0.0120;
                        -0.0017,  0.0003, 0.0003, 0.0003, 0.0000, -0.0010, -0.0038, -0.0082, -0.0140;
                        -0.0007,  0.0003, 0.0003, 0.0003, 0.0000, -0.0017, -0.0048, -0.0099, -0.0161;
                         0.0000,  0.0007, 0.0007, 0.0007, 0.0000, -0.0024, -0.0062, -0.0116, -0.0178;
                         0.0007,  0.0010, 0.0010, 0.0010, 0.0000, -0.0027, -0.0068, -0.0130, -0.0195;
                         0.0017,  0.0017, 0.0017, 0.0017, 0.0000, -0.0031, -0.0082, -0.0144, -0.0209;
                         0.0024,  0.0024, 0.0024, 0.0021, 0.0000, -0.0038, -0.0092, -0.0154, -0.0219;
                         0.0031,  0.0031, 0.0031, 0.0024, 0.0000, -0.0041, -0.0103, -0.0168, -0.0226;
                         0.0041,  0.0041, 0.0041, 0.0027, 0.0000, -0.0048, -0.0110, -0.0175, -0.0229;
                         0.0048,  0.0048, 0.0048, 0.0034, 0.0000, -0.0051, -0.0116, -0.0178, -0.0229];

    % Mach = 10
    Cnel_base(:,:,6) = [-0.0024, -0.0003, 0.0000, 0.0000, 0.0000, -0.0003, -0.0024, -0.0062, -0.0113;
                        -0.0014,  0.0000, 0.0000, 0.0000, 0.0000, -0.0010, -0.0034, -0.0079, -0.0137;
                        -0.0007,  0.0003, 0.0003, 0.0003, 0.0000, -0.0014, -0.0048, -0.0096, -0.0164;
                         0.0000,  0.0007, 0.0007, 0.0007, 0.0000, -0.0021, -0.0058, -0.0116, -0.0192;
                         0.0007,  0.0010, 0.0010, 0.0010, 0.0000, -0.0027, -0.0068, -0.0134, -0.0226;
                         0.0014,  0.0014, 0.0014, 0.0014, 0.0000, -0.0031, -0.0082, -0.0158, -0.0267;
                         0.0021,  0.0021, 0.0021, 0.0017, 0.0000, -0.0038, -0.0096, -0.0185, -0.0305;
                         0.0027,  0.0027, 0.0031, 0.0024, 0.0000, -0.0041, -0.0113, -0.0219, -0.0322;
                         0.0034,  0.0038, 0.0038, 0.0031, 0.0000, -0.0051, -0.0137, -0.0253, -0.0305;
                         0.0041,  0.0048, 0.0051, 0.0034, 0.0000, -0.0062, -0.0161, -0.0267, -0.0281];

    % Mach = 20
    Cnel_base(:,:,7) = [-0.0017, 0.0000, 0.0000, 0.0000, 0.0000,  0.0000, -0.0017, -0.0051, -0.0099;
                        -0.0010, 0.0000, 0.0000, 0.0000, 0.0000, -0.0007, -0.0027, -0.0068, -0.0123;
                        -0.0007, 0.0000, 0.0000, 0.0000, 0.0000, -0.0014, -0.0041, -0.0086, -0.0144;
                         0.0000, 0.0000, 0.0000, 0.0000, 0.0000, -0.0017, -0.0051, -0.0103, -0.0168;
                         0.0007, 0.0007, 0.0007, 0.0007, 0.0000, -0.0021, -0.0062, -0.0120, -0.0182;
                         0.0010, 0.0010, 0.0010, 0.0010, 0.0000, -0.0027, -0.0075, -0.0134, -0.0199;
                         0.0017, 0.0017, 0.0017, 0.0017, 0.0000, -0.0034, -0.0082, -0.0147, -0.0212;
                         0.0024, 0.0024, 0.0024, 0.0021, 0.0000, -0.0041, -0.0096, -0.0158, -0.0219;
                         0.0031, 0.0031, 0.0031, 0.0027, 0.0000, -0.0045, -0.0103, -0.0168, -0.0226;
                         0.0041, 0.0041, 0.0041, 0.0031, 0.0000, -0.0048, -0.0110, -0.0175, -0.0226];
    
    Cnel = interp3(delta_el_base, AoA_base, M_base, Cnel_base, delta_el, AoA, M, 'spline');
end