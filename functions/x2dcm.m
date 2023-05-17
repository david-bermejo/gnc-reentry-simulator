function res = x2dcm(angle)
    % Computes the Direction Cosine Matrix corresponding to a rotation
    % of angle 'angle' (in radians) around the x- axis.
    % 

    sa = sin(angle);
    ca = cos(angle);

    res = [1, 0, 0; 0, ca, sa; 0, -sa, ca];
end