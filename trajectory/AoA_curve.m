function res = AoA_curve(M)
    f = 15;
    g = 40;
    
    K = 1;
    M_c = 7;

    s = (1 + tanh(K.*(M - M_c))) ./ 2;
    res = g.*s + (1-s).*f;
end