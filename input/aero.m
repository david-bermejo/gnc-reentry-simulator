function simdb_aero = aero()
    %% Aerodynamic coefficients database
    % Angle of attack table:
    simdb_aero.AoA_table = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45];          % [deg]
    % Mach number table:
    simdb_aero.M_table = [1.2, 1.5, 2, 3, 5, 10, 20, 30, 40];               % [-]
    % Elevon deflection table:
    simdb_aero.delta_el_table = [-40, -30, -20, -10, 0, 10, 20, 30, 40];    % [deg]
    % Rudder deflection table:
    simdb_aero.delta_rl_table = [0, 10, 20, 30, 40];                        % [deg]
    % Body flap deflection table:
    simdb_aero.delta_b_table = [-20, -10, 0, 10, 20, 30];                   % [deg]
    
    simdb_aero = append2struct(simdb_aero, {'CL0', 'CLb', 'CLel'}, lift_coefficients(), 3);
    simdb_aero = append2struct(simdb_aero, {'CSb0', 'CSbel', 'CSel', 'CSrl'}, sideforce_coefficients(), 4);
    simdb_aero = append2struct(simdb_aero, {'CD0', 'CDh', 'CDrl', 'CDel', 'CDb'}, drag_coefficients(), 5);
    simdb_aero = append2struct(simdb_aero, {'Clb0', 'Clel'}, roll_coefficients(), 2);
    simdb_aero = append2struct(simdb_aero, {'Cm0', 'Cmb', 'Cmel'}, pitch_coefficients(), 3);
    simdb_aero = append2struct(simdb_aero, {'Cnb0', 'Cnbrl', 'Cnel', 'Cnrl'}, yaw_coefficients(), 4);
    simdb_aero = append2struct(simdb_aero, {'dCD_dAoA', 'dCL_dAoA', 'dCm_dAoA', 'dCm_del', 'dCl_del', 'dCn_del', 'dCn_drl'}, coeff_derivatives(simdb_aero), 7);
end