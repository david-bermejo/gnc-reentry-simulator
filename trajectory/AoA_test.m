N = 400;
Ms = linspace(2,20,N);

figure;
plot(Ms, AoA_curve(Ms), 'b');
title('Angle of Attack Profile');
xlabel('Mach number [-]');
ylabel('AoA [deg]');
grid;