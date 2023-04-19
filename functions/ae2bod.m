function res = ae2bod(alpha, beta)
    sa = sin(alpha);
    ca = cos(alpha);
    sb = sin(beta);
    cb = cos(beta);

    res = [ca*cb, -ca*sb, -sa;
              sb,     cb,   0;
           sa*cb, -sa*sb,  ca];
end