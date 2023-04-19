% Cdb = f(AoA, delta_b, M)

function Cdb = drag_inc_body_flap_interpolant(method)
    AoA_base = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]; % deg
    delta_b_base = [-20, -10, 0, 10, 20, 30]; % deg
    M_base = [1.2, 1.5, 2, 3, 5, 10, 20, 40, 60];
    
    % Mach = 1.2
    Cdb_base = [ 0.040,  0.016, 0.000, -0.005, -0.003, 0.005;
                 0.030,  0.010, 0.000,  0.000,  0.004, 0.017;
                 0.018,  0.004, 0.000,  0.016,  0.002, 0.028;
                 0.008,  0.000, 0.000,  0.005,  0.018, 0.040;
                 0.000, -0.003, 0.000,  0.009,  0.027, 0.053;
                -0.005, -0.005, 0.000,  0.013,  0.037, 0.067;
                -0.011, -0.009, 0.000,  0.018,  0.046, 0.078;
                -0.018, -0.013, 0.000,  0.022,  0.053, 0.089;
                -0.026, -0.017, 0.000,  0.026,  0.060, 0.096;
                -0.033, -0.021, 0.000,  0.029,  0.065, 0.100];
    
    % Mach = 1.5
    Cdb_base(:,:,2) = [ 0.021,  0.009, 0.000, -0.003, 0.000, 0.008;
                        0.016,  0.006, 0.000,  0.000, 0.005, 0.017;
                        0.010,  0.002, 0.000,  0.002, 0.011, 0.028;
                        0.004, -0.001, 0.000,  0.005, 0.018, 0.040;
                       -0.001, -0.003, 0.000,  0.009, 0.027, 0.053;
                       -0.006, -0.005, 0.000,  0.013, 0.036, 0.066;
                       -0.011, -0.009, 0.000,  0.018, 0.045, 0.078;
                       -0.018, -0.013, 0.000,  0.022, 0.053, 0.089;
                       -0.025, -0.017, 0.000,  0.026, 0.060, 0.096;
                       -0.033, -0.021, 0.000,  0.029, 0.065, 0.100];
    
    % Mach = 2
    Cdb_base(:,:,3) = [ 0.011,  0.005, 0.000, -0.002, 0.000, 0.008;
                        0.009,  0.003, 0.000,  0.000, 0.005, 0.017;
                        0.006,  0.001, 0.000,  0.002, 0.011, 0.028;
                        0.002,  0.000, 0.000,  0.005, 0.018, 0.040;
                       -0.002, -0.002, 0.000,  0.009, 0.027, 0.053;
                       -0.006, -0.005, 0.000,  0.013, 0.036, 0.066;
                       -0.011, -0.009, 0.000,  0.017, 0.045, 0.078;
                       -0.018, -0.013, 0.000,  0.022, 0.053, 0.089;
                       -0.025, -0.017, 0.000,  0.026, 0.060, 0.096;
                       -0.033, -0.021, 0.000,  0.030, 0.065, 0.100];
    
    % Mach = 3
    Cdb_base(:,:,4) = [ 0.006,  0.003, 0.000, -0.001, 0.001, 0.009;
                        0.004,  0.002, 0.000,  0.000, 0.005, 0.017;
                        0.003,  0.000, 0.000,  0.002, 0.011, 0.028;
                        0.001, -0.001, 0.000,  0.005, 0.018, 0.040;
                       -0.002, -0.002, 0.000,  0.009, 0.027, 0.053;
                       -0.006, -0.005, 0.000,  0.013, 0.036, 0.066;
                       -0.011, -0.009, 0.000,  0.017, 0.045, 0.078;
                       -0.018, -0.013, 0.000,  0.022, 0.053, 0.089;
                       -0.025, -0.017, 0.000,  0.026, 0.060, 0.096;
                       -0.033, -0.021, 0.000,  0.030, 0.065, 0.100];
    
    % Mach = 5
    Cdb_base(:,:,5) = [ 0.002,  0.001, 0.000, 0.000, 0.002, 0.010;
                        0.001,  0.000, 0.000, 0.000, 0.005, 0.017;
                        0.001,  0.000, 0.000, 0.002, 0.011, 0.028;
                        0.000, -0.001, 0.000, 0.005, 0.018, 0.040;
                       -0.002, -0.002, 0.000, 0.009, 0.027, 0.053;
                       -0.006, -0.005, 0.000, 0.013, 0.036, 0.066;
                       -0.011, -0.009, 0.000, 0.017, 0.045, 0.078;
                       -0.018, -0.013, 0.000, 0.022, 0.053, 0.089;
                       -0.025, -0.017, 0.000, 0.026, 0.060, 0.096;
                       -0.033, -0.021, 0.000, 0.030, 0.065, 0.100];

    % Mach = 10
    Cdb_base(:,:,6) = [ 0.000,  0.000, 0.000, 0.000, 0.002, 0.010;
                        0.000,  0.000, 0.000, 0.000, 0.005, 0.017;
                        0.000,  0.000, 0.000, 0.002, 0.011, 0.028;
                       -0.001, -0.001, 0.000, 0.005, 0.018, 0.040;
                       -0.003, -0.003, 0.000, 0.009, 0.027, 0.053;
                       -0.006, -0.005, 0.000, 0.013, 0.036, 0.066;
                       -0.011, -0.009, 0.000, 0.017, 0.045, 0.078;
                       -0.018, -0.013, 0.000, 0.022, 0.053, 0.089;
                       -0.025, -0.017, 0.000, 0.026, 0.060, 0.096;
                       -0.033, -0.021, 0.000, 0.030, 0.065, 0.100];
    
    % Mach = 20
    Cdb_base(:,:,7) = [ 0.000,  0.000, 0.000, 0.000, 0.002, 0.008;
                        0.000,  0.000, 0.000, 0.000, 0.004, 0.015;
                        0.000,  0.000, 0.000, 0.002, 0.009, 0.024;
                       -0.001, -0.001, 0.000, 0.004, 0.015, 0.036;
                       -0.002, -0.002, 0.000, 0.008, 0.024, 0.049;
                       -0.004, -0.004, 0.000, 0.012, 0.033, 0.062;
                       -0.009, -0.007, 0.000, 0.016, 0.043, 0.076;
                       -0.016, -0.011, 0.000, 0.021, 0.051, 0.086;
                       -0.023, -0.016, 0.000, 0.025, 0.059, 0.095;
                       -0.031, -0.020, 0.000, 0.029, 0.065, 0.100];

    % Mach = 30
    Cdb_base(:,:,8) = [ 0.000,  0.000, 0.000, 0.000, 0.002, 0.008;
                        0.000,  0.000, 0.000, 0.000, 0.004, 0.015;
                        0.000,  0.000, 0.000, 0.002, 0.009, 0.024;
                       -0.001, -0.001, 0.000, 0.004, 0.015, 0.036;
                       -0.002, -0.002, 0.000, 0.008, 0.024, 0.049;
                       -0.004, -0.004, 0.000, 0.012, 0.033, 0.062;
                       -0.009, -0.007, 0.000, 0.016, 0.043, 0.076;
                       -0.016, -0.011, 0.000, 0.021, 0.051, 0.086;
                       -0.023, -0.016, 0.000, 0.025, 0.059, 0.095;
                       -0.031, -0.020, 0.000, 0.029, 0.065, 0.100];

    % Mach = 40
    Cdb_base(:,:,9) = [ 0.000,  0.000, 0.000, 0.000, 0.002, 0.008;
                        0.000,  0.000, 0.000, 0.000, 0.004, 0.015;
                        0.000,  0.000, 0.000, 0.002, 0.009, 0.024;
                       -0.001, -0.001, 0.000, 0.004, 0.015, 0.036;
                       -0.002, -0.002, 0.000, 0.008, 0.024, 0.049;
                       -0.004, -0.004, 0.000, 0.012, 0.033, 0.062;
                       -0.009, -0.007, 0.000, 0.016, 0.043, 0.076;
                       -0.016, -0.011, 0.000, 0.021, 0.051, 0.086;
                       -0.023, -0.016, 0.000, 0.025, 0.059, 0.095;
                       -0.031, -0.020, 0.000, 0.029, 0.065, 0.100];
    
    [X, Y, Z] = ndgrid(AoA_base, delta_b_base, M_base);
    Cdb = griddedInterpolant(X, Y, Z, Cdb_base, method);
end