function g = gravity(R, delta)
    J2 = 1960.45e-6; % Coefficient
    Re = 3396.2; % km
    mu = 0.042828e6; % km^3/s^2
    
    g = [mu/R^2*(1 - 1.5*J2*(Re/R)^2*(3*sin(delta)^2 - 1));
         0;
         -3*J2*mu/R^2*(Re/R)^2*sin(delta)*cos(delta)];
end