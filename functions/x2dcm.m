function res = x2dcm(angle)
    sa = sin(angle);
    ca = cos(angle);

    res = [1, 0, 0; 0, ca, sa; 0, -sa, ca];
end