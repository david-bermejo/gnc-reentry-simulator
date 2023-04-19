function [T, p, rho] = atmosphere_interpolant(method)
    h_base = [0.00000e+00, 8.82353e+03, 1.76471e+04, 2.64706e+04, ...
              3.52941e+04, 4.41176e+04, 5.29412e+04, 6.17647e+04, ...
              7.05882e+04, 7.94118e+04, 8.82353e+04, 9.70588e+04, ...
              1.05882e+05, 1.14706e+05, 1.23529e+05, 1.32353e+05, ...
              1.41176e+05, 1.50000e+05, 1.58824e+05, 1.67647e+05, ...
              1.76471e+05, 1.85294e+05, 1.94118e+05, 2.02941e+05, ...
              2.11765e+05, 2.20588e+05, 2.29412e+05, 2.38235e+05, ...
              2.47059e+05, 2.55882e+05, 2.64706e+05, 2.73529e+05, ...
              2.82353e+05, 2.91176e+05, 3.00000e+05];

    T_base = [2.89199e+02, 2.20857e+02, 2.07267e+02, 1.91026e+02, ...
              1.73881e+02, 1.66118e+02, 1.51538e+02, 1.43525e+02, ...
              1.59130e+02, 1.52043e+02, 1.38857e+02, 1.32260e+02, ...
              1.32193e+02, 1.21148e+02, 1.16629e+02, 1.50678e+02, ...
              1.79820e+02, 1.95527e+02, 2.02929e+02, 2.06235e+02, ...
              2.07682e+02, 2.08320e+02, 2.08610e+02, 2.08749e+02, ...
              2.08819e+02, 2.08851e+02, 2.08866e+02, 2.08877e+02, ...
              2.08882e+02, 2.08887e+02, 2.08891e+02, 2.08895e+02, ...
              2.08897e+02, 2.08897e+02, 2.08897e+02];

    p_base = [6.29418e+02, 2.96245e+02, 1.33525e+02, 5.69037e+01, ...
              2.24480e+01, 8.32577e+00, 2.91799e+00, 9.24384e-01, ...
              3.07672e-01, 1.06770e-01, 3.43631e-02, 1.02419e-02, ...
              2.99983e-03, 8.46391e-04, 2.13781e-04, 6.46091e-05, ...
              2.51289e-05, 1.10899e-05, 5.23041e-06, 2.57289e-06, ...
              1.31152e-06, 6.93920e-07, 3.82760e-07, 2.20743e-07, ...
              1.32748e-07, 8.33574e-08, 5.54923e-08, 3.75294e-08, ...
              2.69693e-08, 1.93806e-08, 1.46973e-08, 1.11763e-08, ...
              8.44116e-09, 6.55613e-09, 5.18338e-09];

    rho_base = [1.33165e-02, 7.01575e-03, 3.37000e-03, 1.55890e-03, ...
                6.75779e-04, 2.62303e-04, 1.00871e-04, 3.39167e-05, ...
                1.01321e-05, 3.67188e-06, 1.29206e-06, 4.05599e-07, ...
                1.18539e-07, 3.62459e-08, 9.58984e-09, 2.26016e-09, ...
                7.22396e-10, 2.86819e-10, 1.27684e-10, 6.02866e-11, ...
                2.95394e-11, 1.49383e-11, 7.80681e-12, 4.22195e-12, ...
                2.35004e-12, 1.35369e-12, 8.34978e-13, 5.10420e-13, ...
                3.43851e-13, 2.24152e-13, 1.60917e-13, 1.13939e-13, ...
                7.95688e-14, 5.77415e-14, 4.26350e-14];

    T = griddedInterpolant(h_base, T_base, method);
    p = griddedInterpolant(h_base, p_base, method);
    rho = griddedInterpolant(h_base, rho_base, method);
end