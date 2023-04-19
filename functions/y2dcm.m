function res = y2dcm(angle)
    sa = sin(angle);
    ca = cos(angle);

    res = [ca, 0, -sa; 0, 1, 0; sa, 0, ca];
end