function [CL, CS, CD, Cl, Cm, Cn] = aerodynamic_coefficients(AoA, M, AoS, delta_rl, delta_rr, delta_el, delta_er, delta_b)
    CL = lift_clean(AoA, M) + ...
         lift_inc_left_elevon(AoA, delta_el, M) + ...
         lift_inc_left_elevon(AoA, delta_er, M) + ...
         lift_inc_body_flap(AoA, delta_b, M);
    
    CS = side_inc_left_rudder(AoA, delta_rl, M) - ...
         side_inc_left_rudder(AoA, delta_rr, M) + ...
         side_inc_left_elevon(AoA, delta_el, M) - ...
         side_inc_left_elevon(AoA, delta_el, M) + (...
         side_clean(AoA, M) + ...
         side_inc_deriv_left_elevon(AoA, delta_el, M) + ...
         side_inc_deriv_left_elevon(AoA, delta_er, M) ...
         ) * AoS;
    
    CD = drag_clean(AoA, M) + ...
         drag_inc_left_rudder(AoA, delta_rl, M) + ...
         drag_inc_left_rudder(AoA, delta_rr, M) + ...
         drag_inc_left_elevon(AoA, delta_el, M) + ...
         drag_inc_left_elevon(AoA, delta_er, M) + ...
         drag_inc_body_flap(AoA, delta_b, M);

    Cl = roll_clean(AoA, M) * AoS + ...
         roll_inc_left_elevon(AoA, delta_el, M) - ...
         roll_inc_left_elevon(AoA, delta_er, M);

    Cm = pitch_clean(AoA, M) + ...
         pitch_inc_left_elevon(AoA, delta_el, M) + ...
         pitch_inc_left_elevon(AoA, delta_er, M) + ...
         pitch_inc_body_flap(AoA, delta_b, M);

    Cn = yaw_inc_left_rudder(AoA, delta_rl, M) - ...
         yaw_inc_left_rudder(AoA, delta_rr, M) + ...
         yaw_inc_left_elevon(AoA, delta_el, M) - ...
         yaw_inc_left_elevon(AoA, delta_er, M) + (...
         yaw_clean(AoA, M) + ...
         yaw_inc_deriv_left_rudder(AoA, delta_rl, M) + ...
         yaw_inc_deriv_left_rudder(AoA, delta_rr, M) ...
         ) * AoS;
end