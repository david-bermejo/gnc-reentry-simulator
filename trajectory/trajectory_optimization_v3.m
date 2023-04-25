clear
close all
addpath(pwd + "\..\direct_collocation");
addpath(pwd + "\..\fcn\aerodynamic_coefficients")
addpath(pwd + "\..\fcn\atmospheric_model")
addpath(pwd + "\..\fcn\gravity_model")

%% Parameters setup
params.mu = 0.042828e15;            % [m^3/s^2]
params.Rp = 3396.2 * 1e3;           % [m]
params.omega_cb = 7.07765809e-5;    % [rad/s]
params.gamma = 1.2941;              % [-]
params.Rg = 8314.3/44.01;           % [-]
params.Sref = 110;                  % [m^2]
params.Rn = 1.2;                    % [m]
params.m = 26029.0;                 % [kg]
params.g0 = 9.81;                   % [m/s^2]

params.Qdot_max = 450;              % [kW/m^2]
params.q_max = 10000;               % [Pa]
params.n_max = 2;                   % [-]
params.d_chi_max = deg2rad(4);      % [rad]

params.drag_clean = drag_clean_interpolant('linear');
params.lift_clean = lift_clean_interpolant('linear');
params.pitch_clean = pitch_clean_interpolant('linear');
params.pitch_flap = pitch_inc_body_flap_interpolant('linear');
params.drag_flap = drag_inc_body_flap_interpolant('linear');
params.lift_flap = lift_inc_body_flap_interpolant('linear');
[params.T_interp, params.p_interp, params.rho_interp] = atmosphere('linear');

%% Initial conditions
R0 = params.Rp + 120.0e3;           % [m]
tau0 = deg2rad(22.775755);               % [rad]
delta0 = deg2rad(-17.866444);             % [rad]
V0 = 5411.892434;                          % [m/s]
gamma0 = deg2rad(-8.535284);             % [rad]
chi0 = deg2rad(37.845789);               % [rad]
params.x0 = [R0, tau0, delta0, V0, gamma0, chi0]';

%% Terminal conditions (Perseverance landing site)
Rf = params.Rp + 10e3;               % [m]
tauf = deg2rad(77.45);               % [rad]
deltaf = deg2rad(18.44);             % [rad]
Vf = 800.0;                         % [m/s]
gammaf = deg2rad(-2.0);             % [rad]
chif = deg2rad(90.0);               % [rad]
params.xf = [Rf, tauf, deltaf, Vf, gammaf, chif]';

lb = [5000; deg2rad(-15); deg2rad(0); deg2rad(-60); deg2rad(30)];
ub = [7000; deg2rad(-7); deg2rad(70); deg2rad(5); deg2rad(60)];
options = optimoptions("surrogateopt", "Display", "iter", ...
    "UseParallel", true, 'UseVectorized',true, ...
    'MaxFunctionEvaluations', 250000);

if false
    [x, fval] = surrogateopt(@(x) objval(x, params), lb, ub, options);
    
    fprintf("V0: %f m/s\n", x(1));
    fprintf("gamma0: %f deg\n", rad2deg(x(2)));
    fprintf("tau0: %f deg\n", rad2deg(x(3)));
    fprintf("delta0: %f deg\n", rad2deg(x(4)));
    fprintf("chi0: %f deg\n", rad2deg(x(5)));
    
    %% Obtain reentry solution for optimal initial conditions
    % x = [V0, gamma0, tau0, delta0, chi0];
    [ts, y, u, ~] = objfun(x, params);
    
    %% Save variables for Simulink
    
    %e0 = y(4,1)^2 / 2 - params.mu ./ y(1,1);
    %energy = (y(4,:).^2 ./ 2 - params.mu ./ y(1,:)) ./ e0;
    e0 = y(4,1)^2/2 + params.mu/y(1,1)^2*(y(1,1) - params.Rp);
    energy = y(4,:).^2./2 + params.mu./y(1,:).^2.*(y(1,:) - params.Rp);
    energy = energy / e0;
    
    tau0 = x(3);
    delta0 = x(4);
    V0 = x(1);
    gamma0 = x(2);
    chi0 = x(5);
    
    save('../input/trajectory.mat', 'ts', 'energy', 'y', 'u', ...
        'R0', 'tau0', 'delta0', 'V0', 'gamma0', 'chi0');
else
    load('../input/trajectory.mat', 'ts', 'energy', 'y', 'u', ...
        'R0', 'tau0', 'delta0', 'V0', 'gamma0', 'chi0');
end

%% Plot results
figure();
plot(ts, energy);
title('Energy');
xlabel('t [s]');
ylabel('Energy');
grid;

figure();
plot(energy, rad2deg(u(1,:)));
title('Control Law: u(t)');
xlabel('Energy [-]');
ylabel('Bank Angle [deg]');
grid;

figure();
plot(ts, rad2deg(u(1,:)));
title('Control Law: u(t)');
xlabel('t [s]');
ylabel('Bank Angle [deg]');
grid;

figure();
plot(ts, u(2,:));
title('Body flap: \delta_b(t)');
xlabel('t [s]');
ylabel('Body flap [deg]');
grid;

figure;
plot(ts, (y(1,:) - params.Rp) ./ 1000);
title('Reentry Envelope');
xlabel('t [s]');
ylabel('Height [km]');
grid;

figure;
plot(ts, y(4,:) ./ 1000);
title('Velocity');
xlabel('t [s]');
ylabel('Velocity [km/s]');
grid;

figure;
plot(y(4,:) ./ 1000, (y(1,:) - params.Rp) ./ 1000);
title('Flight Envelope');
xlabel('Velocity [km/s]');
ylabel('Height [km]');
grid;

figure;
subplot(2, 1, 1);
plot(ts, rad2deg(y(2,:)));
title('Longitude');
xlabel('t [s]');
ylabel('\tau [deg]');
grid;
subplot(2, 1, 2);
plot(ts, rad2deg(y(3,:)));
title('Latitude');
xlabel('t [s]');
ylabel('\delta [deg]');
grid;

figure;
subplot(2, 1, 1);
plot(ts, rad2deg(y(5,:)));
title('Flight Path Angle');
xlabel('t [s]');
ylabel('\gamma [deg]');
grid;
subplot(2, 1, 2);
plot(ts, rad2deg(y(6,:)));
hold on;
plot(ts, rad2deg(delta_heading(y, params) + y(6,:)));
plot(ts, rad2deg(delta_heading(y, params)));
hold off;
title('Heading');
xlabel('t [s]');
ylabel('\chi [deg]');
legend('\chi', '\chi_t', '\Delta\chi');
grid;

figure;
plot(ts, rad2deg(delta_heading(y, params)));
title('Heading angle error \Delta\chi(t)');
xlabel('t [s]');
ylabel('\Delta\chi [deg]');
grid;

constr = calc_constraints(ts, y, u, params);
figure;
plot(ts, constr(3,:));
title('Normal Acceleration');
xlabel('t [s]');
ylabel('n_{max} [-]');
grid;

figure;
plot(ts, constr(1,:));
title('Heat Flux profile');
xlabel('t [s]');
ylabel('q_{tot} [kW/m^2]');
grid;

h = linspace(R0-params.Rp, Rf-params.Rp, 100);
V1 = zeros(1, 100);
V2 = zeros(1, 100);
V3 = zeros(1, 100);
V4 = zeros(1, 100);


for i=1:100
    V1(i) = fzero(@(x) glide_eqn(40.0, x, h(i), params), 3000);
    V2(i) = fzero(@(x) total_heat_rate_eqn(x, h(i), params), 3000);
    V3(i) = g_load_eqn(40.0, h(i), params);
    V4(i) = qdyn_eqn(h(i), params);
end

min_data = min(V2, V3);
x2 = [min_data, fliplr(V1)]/1000;
inBetween = [h, fliplr(h)]/1000;

figure;
fill(x2, inBetween, 'y');
hold on;
plot(y(4,:)/1000, (y(1,:) - params.Rp)/1000, 'LineWidth', 1.5);
plot(V1/1000, h/1000);
plot(V2/1000, h/1000, '--');
plot(V3/1000, h/1000, '--');
plot(V4/1000, h/1000, '--');
hold off;
xlim([0, V0/1000*1.025]);
title('Entry Corridor');
xlabel('V [km/s]');
ylabel('h [km]');
name = legend('', 'Numerical Solution', 'Glide, $\alpha = 40^{\circ}$', ...
        "$\dot{Q}_{max} = " + params.Qdot_max + "\;kW/m^2$", ...
        "$n_{max} = " + params.n_max + "$, $\alpha = 40^{\circ}$", ...
        "$q_{max} = " + params.q_max + "\;Pa$");
set(name,'Interpreter','latex');
grid;

figure;
plot(rad2deg(y(2,:)), rad2deg(y(3,:)), 'y', 'LineWidth',1.5);
hold on
plot(rad2deg(tauf), rad2deg(deltaf), 'pentagram', 'Color','white', 'MarkerSize',7, 'LineWidth',2);
axis tight;
[img, cmap] = imread('map.png', 'png');
h = image([-180, 180], -[-90, 90], img);
uistack(h, "bottom");
colormap(cmap);
pbaspect([2 1 1]);
title('Reentry Ground Track');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
xticks(-180:30:180);
xtickformat('degrees');
yticks(-90:15:90);
ytickformat('degrees');

%% Equations
function J = objval(opt, p)
    [~, ~, ~, J] = objfun(opt, p);
end

function [ts, y, u, J] = objfun(opt, p)
    p.x0(2) = opt(3);
    p.x0(3) = opt(4);
    p.x0(4) = opt(1);
    p.x0(5) = opt(2);
    p.x0(6) = opt(5);

    tspan = [0, 2000];
    ts = tspan(1):1:tspan(end);
    N = length(ts);
    
    u = zeros(2,N);
    y = zeros(6,N);
    y(:,1) = p.x0;

    options = odeset('MaxStep',2, 'Events', @(t, x) myEvent(t, x, p));
    
    for i=2:N
        h = y(1,i-1) - p.Rp;
    
        if h <= 120e3
            T = p.T_interp(h);
            rho = p.rho_interp(h);
            g = gravity(y(1,i-1));
        
            M = y(4,i-1) ./ sqrt(p.gamma.*p.Rg.*T);
            AoA = AoA_curve(M);
        
            Cm0 = p.pitch_clean(AoA, M);
            u(2,i) = fzero(@(x) Cm0 + p.pitch_flap(AoA,x,M), u(2,i-1));

            q_inf = 0.5.*rho.*y(4,i-1).^2;
            CL0 = p.lift_clean(AoA, M);
            CLb = p.lift_flap(AoA, u(2,i), M);
            CL = CL0 + CLb;
            L = q_inf.*CL.*p.Sref;
            
            var = (g - y(4,i-1).^2./y(1,i-1)) .* cos(y(5,i-1)) .* p.m ./ L;
            tmp = min(max(var, -1.0), 1.0);
            
            d_chi = delta_heading(y(:,i-1), p);
            u(1,i) = acos(tmp) .* sgn(d_chi, p.d_chi_max, u(1,i-1));
            u(1,i) = max(min(u(1,i), deg2rad(165)), deg2rad(-165));
        end
    
        % Solve the problem
        solution = ode45(@(t, x) dynamics_fixed(t, x, u(:,i), p), [ts(i-1), ts(i)], y(:,i-1), options);
    
        % Save solution
        y(:,i) = solution.y(:,end);
    
        if y(1,i) <= p.xf(1)
            ts = ts(1:i);
            y = y(:,1:i);
            u = u(:,1:i);
            break;
        end
    end

    if ts(end) == tspan(end)
        J = Inf;
        return;
    end

    [Cin, ~] = nlcon(ts, y, u, p);
    cst = pathCst(Cin);
    cost = objBnd(y(:,1), ts(1), y(:,end), ts(end), p);

    J = cost + 1e6 .* sum((ts(2:end) - ts(1:end-1))./2 .* (cst(:,1:end-1) + cst(:,2:end)), "all");
end

function cost = objBnd(x0, t0, xf, tf, p)
    cost = 10*((xf(2) - p.xf(2))^2 + (xf(3) - p.xf(3))^2);
end

function res = pathCst(x)
    res = (x > 0) .* x.^2;
end

function [Cin, Ceq] = nlcon(t, x, u, p)
    h = x(1,:) - p.Rp;
    T = p.T_interp(h);
    rho = p.rho_interp(h);
    M = x(4,:) ./ sqrt(p.gamma.*p.Rg.*T);
    AoA = AoA_curve(M);

    q_inf = 0.5.*rho.*x(4,:).^2;
    CD0 = p.drag_clean(AoA, M);
    CL0 = p.lift_clean(AoA, M);
    CDb = p.drag_flap(AoA, u(2,:), M);
    CLb = p.lift_flap(AoA, u(2,:), M);
    CD = CD0 + CDb;
    CL = CL0 + CLb;
    D = q_inf.*CD.*p.Sref;
    L = q_inf.*CL.*p.Sref;
    
    Qdot = 1.9027e-7 .* sqrt(rho./p.Rn) .* x(4,:).^3;   % [kW/m^2]
    %AoA_rad = deg2rad(AoA);
    n = sqrt(D.^2 + L.^2) ./ (p.m.*p.g0);     % [-]

    Cin = [Qdot - p.Qdot_max; q_inf - p.q_max; n - p.n_max];
    Ceq = [];
    return 
end

function C = calc_constraints(t, x, u, p)
    [Cin, ~] = nlcon(t, x, u, p);
    C = Cin + [p.Qdot_max; p.q_max; p.n_max];
end

function [value, isterminal, direction] = myEvent(t, x, p)
    value      = x(1) - p.xf(1);
    isterminal = 1;   % Stop the integration
    direction  = -1;
end