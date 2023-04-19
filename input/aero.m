function out = aero()
    %% Aerodynamic coefficients - database initialization
    out.AoA_table = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45];          % deg
    out.M_table = [1.2, 1.5, 2, 3, 5, 10, 20, 30, 40];                       % deg
    out.delta_el_table = [-40, -30, -20, -10, 0, 10, 20, 30, 40];    % deg
    out.delta_rl_table = [0, 10, 20, 30, 40];                        % deg
    out.delta_b_table = [-20, -10, 0, 10, 20, 30];                   % deg
    
    out = append2struct(out, {'CL0', 'CLb', 'CLel'}, lift_coefficients(), 3);
    out = append2struct(out, {'CSb0', 'CSbel', 'CSel', 'CSrl'}, sideforce_coefficients(), 4);
    out = append2struct(out, {'CD0', 'CDh', 'CDrl', 'CDel', 'CDb'}, drag_coefficients(), 5);
    out = append2struct(out, {'Clb0', 'Clel'}, roll_coefficients(), 2);
    out = append2struct(out, {'Cm0', 'Cmb', 'Cmel'}, pitch_coefficients(), 3);
    out = append2struct(out, {'Cnb0', 'Cnbrl', 'Cnel', 'Cnrl'}, yaw_coefficients(), 4);
end