clear
close all
addpath(pwd + "\..\fcn\aerodynamic_coefficients")
addpath(pwd + "\..\fcn\atmospheric_model")

%% Parameters setup
params.mu        = 0.042828e15;     % [m^3/s^2]
params.Rp        = 3396.2 * 1e3;    % [m]
params.omega_cb  = 7.07765809e-5;   % [rad/s]
params.gamma     = 1.2941;          % [-]
params.Rg        = 8314.3/44.01;    % [-]
params.Sref      = 110;             % [m^2]
params.Rn        = 1.2;             % [m]
params.m         = 26029.0;         % [kg]
params.g0        = 9.81;            % [m/s^2]

params.Kq        = 1.9027e-7;       % [kW/m^2]
params.R_max     = params.Rp+150e3; % [m]
params.Qdot_max  = 500;             % [kW/m^2]
params.q_max     = 10000;           % [Pa]
params.n_max     = 2.5;             % [-]
params.gamma_max = deg2rad(0);      % [rad]

params.d_chi0    = deg2rad(4);      % [rad]
params.d_chif    = deg2rad(2.5);    % [rad]
params.M0        = 40;              % [-]
params.Mf        = 5;               % [-]
params.sigma_max = deg2rad(135);    % [rad]
params.sigma_dot = deg2rad(15);     % [rad/s]
params.sigma_dmp = 1.0;             % [-]
params.r         = 0.95;            % [-]

params.drag_clean = drag_clean_interpolant('linear');
params.lift_clean = lift_clean_interpolant('linear');
params.pitch_clean = pitch_clean_interpolant('linear');
params.pitch_flap = pitch_inc_body_flap_interpolant('linear');
params.drag_flap = drag_inc_body_flap_interpolant('linear');
params.lift_flap = lift_inc_body_flap_interpolant('linear');
[params.T_interp, params.p_interp, params.rho_interp] = atmosphere_interpolant('linear');

N = 1000;
h = linspace(0, 150e3, N);
rho = params.rho_interp(h);
T = params.T_interp(h);

figure;
subplot(2,1,1);
semilogx(rho, h/1000, 'b','LineWidth',1.5);
xlabel('Density [kg/m^3]');
ylabel('Height [km]')
grid;
subplot(2,1,2);
plot(T, h/1000, 'b','LineWidth',1.5);
xlabel('Temperature [K]');
ylabel('Height [km]')
grid;

%% Initial conditions
R0 = params.Rp + 120.0e3;           % [m]
tau0 = deg2rad(22.775755);          % [rad]
delta0 = deg2rad(-17.866444);       % [rad]
V0 = 5411.892434;                   % [m/s]
gamma0 = deg2rad(-8.535284);        % [rad]
chi0 = deg2rad(37.845789);          % [rad]
params.x0 = [R0, tau0, delta0, V0, gamma0, chi0]';

%% Terminal conditions (Perseverance landing site)
Rf = params.Rp + 15e3;              % [m]
tauf = deg2rad(77.45);              % [rad]
deltaf = deg2rad(18.44);            % [rad]
Vf = 800.0;                         % [m/s]
gammaf = deg2rad(-2.0);             % [rad]
chif = deg2rad(90.0);               % [rad]
params.xf = [Rf, tauf, deltaf, Vf, gammaf, chif]';

lb = [5600; deg2rad(-14); deg2rad(10); deg2rad(-15); deg2rad(50)];
ub = [6000; deg2rad(-9); deg2rad(60); deg2rad(20); deg2rad(90)];
options = optimoptions("surrogateopt",...
    'PlotFcn','surrogateoptplot',...
    "Display", "iter",...
    "UseParallel", true,...
    'UseVectorized',true,...
    'MaxFunctionEvaluations',250000);

if true
    [x, fval] = surrogateopt(@(x) objval(x, params), lb, ub, options);
    res = objval(x, params);

    fprintf('\nObjective function value; %g\n', res.Fval);
    fprintf('Inequalities: [%f, %f, %f, %f]\n\n', res.Ineq);
    fprintf("V0: %f m/s\n", x(1));
    fprintf("gamma0: %f deg\n", rad2deg(x(2)));
    fprintf("tau0: %f deg\n", rad2deg(x(3)));
    fprintf("delta0: %f deg\n", rad2deg(x(4)));
    fprintf("chi0: %f deg\n", rad2deg(x(5)));
    
    %% Obtain reentry solution for optimal initial conditions
    % x = [V0, gamma0, tau0, delta0, chi0];
    [ts, y, u, bank] = objfun(x, params);
    
    %% Save variables for Simulink
    
    e0 = y(4,1)^2/2 + params.mu/y(1,1)^2*(y(1,1) - params.Rp);
    ef = y(4,end)^2/2 + params.mu/y(1,end)^2*(y(1,end) - params.Rp);
    energy = y(4,:).^2./2 + params.mu./y(1,:).^2.*(y(1,:) - params.Rp);

    tau0 = x(3);
    delta0 = x(4);
    V0 = x(1);
    gamma0 = x(2);
    chi0 = x(5);
    
    save('trajectory.mat', 'ts', 'energy', 'y', 'u',...
        'e0', 'ef', 'R0', 'tau0', 'delta0', 'V0', 'gamma0', 'chi0');
else
    load('trajectory-v33.mat', 'ts', 'energy', 'y', 'u',...
        'R0', 'tau0', 'delta0', 'V0', 'gamma0', 'chi0');
    x = [V0, gamma0, tau0, delta0, chi0];

    [~, ~, ~, bank] = objfun(x, params);
    res = objval(x, params);

    fprintf('\nObjective function value; %g\n', res.Fval);
    fprintf('Inequalities: [%f, %f, %f, %f]\n\n', res.Ineq);

    fprintf("h0: %f km\n", (params.x0(1)-params.Rp)/1000);
    fprintf("tau0: %f deg\n", rad2deg(x(3)));
    fprintf("delta0: %f deg\n", rad2deg(x(4)));
    fprintf("V0: %f m/s\n", x(1));
    fprintf("gamma0: %f deg\n", rad2deg(x(2)));
    fprintf("chi0: %f deg\n\n", rad2deg(x(5)));

    fprintf("hf: %f km\n", (y(1,end)-params.Rp)/1000);
    fprintf("tauf: %f deg\n", rad2deg(y(2,end)));
    fprintf("deltaf: %f deg\n", rad2deg(y(3,end)));
    fprintf("Vf: %f m/s\n", y(4,end));
    fprintf("gammaf: %f deg\n", rad2deg(y(5,end)));
    fprintf("chif: %f deg\n\n", rad2deg(y(6,end)));
end

T = params.T_interp(y(1,:) - params.Rp);
M = y(4,:) ./ sqrt(params.gamma .* params.Rg .* T);
AoA_cmd = AoA_curve(M);

%% Plot results
figure();
subplot(2,1,1);
plot(M, AoA_cmd, 'b','LineWidth',1.5);
title('Angle of Attack profile');
xlabel('M [-]');
ylabel('\alpha [deg]');
grid;
subplot(2,1,2);
plot(ts, rad2deg(bank), 'b','LineWidth',1.5);
title('Bank angle profile');
xlabel('t [s]');
ylabel('\sigma [deg]');
ylim([-150,150]);
grid;

figure();
plot(ts, energy);
title('Energy');
xlabel('t [s]');
ylabel('Energy');
grid;

figure();
plot(ts, rad2deg(u(1,:)), 'b','LineWidth',1.5);
title('Raw Control Law: u(t)');
xlabel('t [s]');
ylabel('Bank Angle [deg]');
grid;

figure();
plot(ts, rad2deg(bank), 'b','LineWidth',1.5);
title('Control Law: u(t)');
xlabel('t [s]');
ylabel('Bank Angle [deg]');
grid;

figure();
plot(ts, u(2,:), 'b','LineWidth',1.5);
title('Body flap: \delta_b(t)');
xlabel('t [s]');
ylabel('Body flap [deg]');
grid;

figure;
plot(ts, (y(1,:) - params.Rp) ./ 1000, 'b','LineWidth',1.5);
title('Reentry Envelope');
xlabel('t [s]');
ylabel('Height [km]');
grid;

figure;
plot(ts, y(4,:) ./ 1000, 'b','LineWidth',1.5);
title('Velocity');
xlabel('t [s]');
ylabel('Velocity [km/s]');
grid;

figure;
plot(y(4,:) ./ 1000, (y(1,:) - params.Rp) ./ 1000, 'b','LineWidth',1.5);
title('Flight Envelope');
xlabel('Velocity [km/s]');
ylabel('Height [km]');
grid;

figure;
subplot(2, 1, 1);
plot(ts, rad2deg(y(2,:)), 'b','LineWidth',1.5);
title('Longitude');
xlabel('t [s]');
ylabel('\tau [deg]');
grid;
subplot(2, 1, 2);
plot(ts, rad2deg(y(3,:)), 'b','LineWidth',1.5);
title('Latitude');
xlabel('t [s]');
ylabel('\delta [deg]');
grid;

figure;
subplot(2, 1, 1);
plot(ts, rad2deg(y(5,:)), 'b','LineWidth',1.5);
title('Flight Path Angle');
xlabel('t [s]');
ylabel('\gamma [deg]');
grid;
subplot(2, 1, 2);
plot(ts, rad2deg(delta_heading(y, params)), 'b','LineWidth',1.5);
title('Heading Angle Error');
xlabel('t [s]');
ylabel('\Delta\chi [deg]');
ylim([-8, 8]);
grid;

figure;
plot(ts, rad2deg(y(6,:)), 'b','LineWidth',1.5);
hold on;
plot(ts, rad2deg(delta_heading(y, params) + y(6,:)), 'LineWidth',1.5);
plot(ts, rad2deg(delta_heading(y, params)), 'LineWidth',1.5);
hold off;
title('Heading');
xlabel('t [s]');
ylabel('\chi [deg]');
legend('\chi', '\chi_t', '\Delta\chi');
grid;

% Plot constraints
constr = calc_constraints(ts, y, u, params);
figure;
subplot(3,1,1);
plot(ts, constr(3,:), 'b','LineWidth',1.5);
title('Normal Acceleration');
xlabel('t [s]');
ylabel('n [-]');
grid;
subplot(3,1,2);
plot(ts, constr(1,:), 'b','LineWidth',1.5);
title('Heat Rate');
xlabel('t [s]');
ylabel('q_t [kW/m^2]');
grid;
subplot(3,1,3);
plot(ts, constr(2,:), 'b','LineWidth',1.5);
title('Dynamic Pressure');
xlabel('t [s]');
ylabel('q_{\infty} [Pa]');
grid;

h = linspace(R0-params.Rp, Rf-params.Rp, 100);
V1 = zeros(1, 100);
V2 = zeros(1, 100);
V3 = zeros(1, 100);
V4 = zeros(1, 100);


for i=1:100
    V1(i) = fzero(@(x) glide_eqn(45.0, x, h(i), params), 3000);
    V2(i) = fzero(@(x) total_heat_rate_eqn(x, h(i), params), 3000);
    V3(i) = fzero(@(x) g_load_eqn(45.0, x, h(i), params), 2000);
    V4(i) = qdyn_eqn(h(i), params);
end

min_data = min(V2, V3);
x2 = [min_data, fliplr(V1)]/1000;
inBetween = [h, fliplr(h)]/1000;

figure;
fill(x2, inBetween, 'y');
hold on;
%plot(y(4,:)/1000, (y(1,:) - params.Rp)/1000, 'b','LineWidth',1.5);
plot(V1/1000, h/1000);
plot(V2/1000, h/1000, '--');
plot(V3/1000, h/1000, '--');
plot(V4/1000, h/1000, '--');
hold off;
xlim([0, V0/1000*1.025]);
title('Entry Corridor');
xlabel('V [km/s]');
ylabel('h [km]');
name = legend('', 'Glide, $\alpha = 45^{\circ}$', ...
        "$\dot{Q}_{max} = " + params.Qdot_max + "\;kW/m^2$", ...
        "$n_{max} = " + params.n_max + "$, $\alpha = 45^{\circ}$", ...
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
yticks(-90:30:90);
ytickformat('degrees');

%% Equations
function f = objval(opt, p)
    [~, ~, ~, ~, Fval, Ineq] = objfun(opt, p);
    f.Fval = Fval;
    f.Ineq = Ineq;
end

function [ts, y, u, bank, Fval, Ineq] = objfun(opt, p)
    p.x0(2) = opt(3);
    p.x0(3) = opt(4);
    p.x0(4) = opt(1);
    p.x0(5) = opt(2);
    p.x0(6) = opt(5);

    N = 5501;
    tspan = [0, 1100];
    ts = linspace(tspan(1), tspan(2), N);
    dt = ts(2) - ts(1);
    
    bank = zeros(1,N);
    u = zeros(2,N);
    y = zeros(6,N);
    y(:,1) = p.x0;

    omega_n = 10;
    damp    = p.sigma_dmp;
    A = [0 1; -omega_n^2 -2*damp*omega_n];
    B = [0; omega_n^2];
    Y = [0; 0];

    options = odeset('MaxStep',2, 'Events', @(t, x) myEvent(t, x, p));
    
    for i=2:N
        [u(1,i), u(2,i)] = bank_angle_cmd(y(:,i-1), u(:,i-1), p);
        if i==2
            Y(1) = u(1,i);
        end

        M = 3;
        for j=1:M
            dY = A*Y + B*u(1,i);
            dY(1) = max(min(dY(1), p.sigma_dot), -p.sigma_dot);
            Y = Y + dY*(dt/M);
        end
        bank(i) = Y(1);
    
        % Solve the problem
        solution = ode45(@(t, x) dynamics(t, x, [bank(i); u(2,i)], p), [ts(i-1), ts(i)], y(:,i-1), options);
    
        % Save solution
        y(:,i) = solution.y(:,end);
    
        if y(1,i) <= p.xf(1)
            ts = ts(1:i);
            y = y(:,1:i);
            u(:,1) = u(:,2);
            u = abs(u(:,1:i));
            bank = bank(1:i);
            bank(1) = bank(2);
            break;

        elseif y(1,i) >= p.R_max
            ts = ts(1:i);
            y = y(:,1:i);
            u(:,1) = u(:,2);
            u = abs(u(:,1:i));
            bank = bank(1:i);
            bank(1) = bank(2);
            Fval = Inf;
            Ineq = ones(4,1)*Inf;
            return;

        elseif i == N
            bank(1) = bank(2);
            u(:,1) = u(:,2);
            u = abs(u);
        end
    end

    [Cin, ~] = nlcon(ts, y, u, p);
    cst = pathCst(Cin);
    cost = objBnd(y(:,1), ts(1), y(:,end), ts(end), p);

    Fval = cost;
    if i < N
        Ineq = sum((ts(2:end) - ts(1:end-1))./2 .* (cst(:,1:end-1) + cst(:,2:end)), 2);
    else
        Ineq = ones(size(cst(:,1)))*Inf;
    end
end

function cost = objBnd(x0, t0, xf, tf, p)
    cost = log10((xf(2) - p.xf(2))^2 + (xf(3) - p.xf(3))^2);
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
    
    Qdot = p.Kq .* sqrt(rho./p.Rn) .* x(4,:).^3;   % [kW/m^2]
    n = sqrt(D.^2 + L.^2) ./ (p.m.*p.g0);     % [-]

    Cin = [Qdot./p.Qdot_max-1; q_inf./p.q_max-1; n./p.n_max-1; x(5,:)-p.gamma_max];
    Ceq = [];
    return 
end

function C = calc_constraints(t, x, u, p)
    [Cin, ~] = nlcon(t, x, u, p);
    C = (Cin + [1;1;1;p.gamma_max]) .* [p.Qdot_max; p.q_max; p.n_max; 1];
end

function [value, isterminal, direction] = myEvent(t, x, p)
    value      = x(1) - p.xf(1);
    isterminal = 1;   % Stop the integration
    direction  = -1;
end
