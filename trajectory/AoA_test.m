N = 100000;
Ms = linspace(2,20,N);

figure;
plot(Ms, AoA_curve(Ms), 'b','LineWidth',1.5);
title('Angle of Attack Profile');
xlabel('Mach number [-]');
ylabel('AoA [deg]');
grid;