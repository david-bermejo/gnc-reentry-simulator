function res = z2dcm(angle)
    % Computes the Direction Cosine Matrix corresponding to a rotation
    % of angle 'angle' (in radians) around the z- axis.
    %

    sa = sin(angle);
    ca = cos(angle);

    res = [ca, sa, 0; -sa, ca, 0; 0, 0, 1];
end