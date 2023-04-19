function res = z2dcm(angle)
    sa = sin(angle);
    ca = cos(angle);

    res = [ca, sa, 0; -sa, ca, 0; 0, 0, 1];
end