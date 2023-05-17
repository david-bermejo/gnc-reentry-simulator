function res = quatconj(q)
    % Computes the conjugate of the given quaternion:
    % It is calculated as follows: q* = [q_w; -q_x; -q_y; -q_z]

    res = [q(1); -q(2:4)];
end