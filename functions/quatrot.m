function res = quatrot(q, v)
    % Rotates the vector v with the provided quaternion q.
    % The quaternion needs to represent a valid rotation (norm(q) == 1)
    
    t = 2*cross(q(2:4), v);
    res = v + q(1)*t + cross(q(2:4), t);
end