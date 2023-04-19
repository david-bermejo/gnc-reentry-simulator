function res = quatrot(q, v)
    t = 2*cross(q(2:4), v);
    res = v + q(1)*t + cross(q(2:4), t);
end