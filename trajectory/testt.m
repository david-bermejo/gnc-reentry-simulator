omega_n = 10;
damp    = 1.0;
lim = deg2rad(10);
A = [0 1; -omega_n^2 -2*damp*omega_n];
B = [0; omega_n^2];
Y = [0; 0];

N = 10000;
ts = linspace(0, 1000, N);
dt = ts(2) - ts(1);
bank = zeros(1,N);

for i=2:N
    if ts(i) > 250 && ts(i) < 500
        u = deg2rad(100);
    elseif ts(i) >= 500
        u = -deg2rad(100);
    else
        u = 0;
    end

    M = 10;
    for j=1:M
        dY = A*Y + B*u;
        dY(1) = max(min(dY(1), lim), -lim);
        Y = Y + dY*(dt/M);
    end
    bank(i) = Y(1);
end

figure
plot(ts, rad2deg(bank))