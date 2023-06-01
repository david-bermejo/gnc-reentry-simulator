function res = AoA_curve(M)
    f = 25;
    g = 45;
    
    K = 2;
    M_c = 9;

    s = (1 + tanh(K.*(M - M_c))) ./ 2;
    res = g.*s + (1-s).*f;
end