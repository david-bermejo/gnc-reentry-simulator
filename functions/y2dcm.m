function res = y2dcm(angle)
    % Computes the Direction Cosine Matrix corresponding to a rotation
    % of angle 'angle' (in radians) around the y- axis.
    %

    sa = sin(angle);
    ca = cos(angle);

    res = [ca, 0, -sa; 0, 1, 0; sa, 0, ca];
end