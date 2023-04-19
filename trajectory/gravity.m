function g = gravity(R)
    mu = 0.042828e15;    % [m^3/s^2]
    g = mu./R.^2;
end