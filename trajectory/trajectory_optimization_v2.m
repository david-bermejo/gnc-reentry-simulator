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
params.q_max = 10000;                % [Pa]
params.n_max = 2;                    % [-]
params.d_chi_max = deg2rad(4);       % [rad]

params.drag_clean = drag_clean_interpolant('spline');
params.lift_clean = lift_clean_interpolant('spline');
params.pitch_clean = pitch_clean_interpolant('spline');
params.pitch_flap = pitch_inc_body_flap_interpolant('spline');
params.drag_flap = drag_inc_body_flap_interpolant('spline');
params.lift_flap = lift_inc_body_flap_interpolant('spline');
[params.T_interp, params.p_interp, params.rho_interp] = atmosphere('spline');

%% Initial conditions
R0 = params.Rp + 120.0e3;           % [m]
tau0 = deg2rad(11.505947);               % [rad]
delta0 = deg2rad(-17.071250);             % [rad]
V0 = 5634.742320;                          % [m/s]
gamma0 = deg2rad(-8.956112);             % [rad]
chi0 = deg2rad(60);               % [rad]
params.x0 = [R0, tau0, delta0, V0, gamma0, chi0]';

%% Terminal conditions
Rf = params.Rp + 10e3;               % [m]
tauf = deg2rad(77.45);               % [rad]
deltaf = deg2rad(18.44);             % [rad]
Vf = 800.0;                         % [m/s]
gammaf = deg2rad(-2.0);             % [rad]
chif = deg2rad(90.0);               % [rad]
params.xf = [Rf, tauf, deltaf, Vf, gammaf, chif]';

%% Calculate the initial guess
tspan = [0, 1500];
ts = tspan(1):1:tspan(end);
N = length(ts);

u = zeros(2,N);
y = zeros(6,N);
y(:,1) = params.x0;

options = odeset('MaxStep',2, 'Events', @(t, x) myEvent(t, x, params));

for i=2:N
    h = y(1,i-1) - params.Rp;

    if h <= 90e3
        T = params.T_interp(h);
        rho = params.rho_interp(h);
        g = gravity(y(1,i-1));
    
        M = y(4,i-1) ./ sqrt(params.gamma.*params.Rg.*T);
        AoA = AoA_curve(M);
    
        Cm0 = params.pitch_clean(AoA, M);
        u(2,i) = fzero(@(x) Cm0 + params.pitch_flap(AoA,x,M), u(2,i-1));

        q_inf = 0.5.*rho.*y(4,i-1).^2;
        CL0 = params.lift_clean(AoA, M);
        CLb = params.lift_flap(AoA, u(2,i), M);
        CL = CL0 + CLb;
        L = q_inf.*CL.*params.Sref;
        
        var = (g - y(4,i-1).^2./y(1,i-1)) .* cos(y(5,i-1)) .* params.m ./ L;
        tmp = min(max(var, -1.0), 1.0);
        
        d_chi = delta_heading(y(:,i-1), params);
        u(1,i) = acos(tmp) .* sgn(d_chi, params.d_chi_max, u(1,i-1));
    end

    % Solve the problem
    solution = ode45(@(t, x) dynamics_fixed(t, x, u(:,i), params), [ts(i-1), ts(i)], y(:,i-1), options);

    % Save solution
    y(:,i) = solution.y(:,end);

    if y(1,i) <= params.xf(1)
        ts = ts(1:i);
        y = y(:,1:i);
        u = u(:,1:i);
        break;
    end
end

figure();
plot(ts, rad2deg(u(1,:)));
title('Control Law: u(t)');
xlabel('t [s]');
ylabel('Bank Angle [deg]');
grid;

figure();
plot(ts, u(2,:));
title('Body flap: $\delta_b$(t)');
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
xlabel('t [s]');
ylabel('long [deg]');
grid;
subplot(2, 1, 2);
plot(ts, rad2deg(y(3,:)));
xlabel('t [s]');
ylabel('lat [deg]');
grid;

figure;
subplot(2, 1, 1);
plot(ts, rad2deg(y(5,:)));
xlabel('t [s]');
ylabel('gamma [deg]');
grid;
subplot(2, 1, 2);
plot(ts, rad2deg(y(6,:)));
hold on;
plot(ts, rad2deg(delta_heading(y, params) + y(6,:)));
plot(ts, rad2deg(delta_heading(y, params)));
hold off;
xlabel('t [s]');
ylabel('chi [deg]');
legend('chi', 'chi_d', 'delta_chi');
grid;

figure;
plot(ts, rad2deg(delta_heading(y, params)));
title('Heading angle error delta_chi(t)');
xlabel('t [s]');
ylabel('Heading Angle Error [deg]');
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
V5 = zeros(1, 100);

for i=1:100
    V1(i) = fzero(@(x) glide_eqn(40.0, x, h(i), params), 3000);
    V2(i) = fzero(@(x) glide_eqn(15.0, x, h(i), params), 3000);
    V3(i) = fzero(@(x) total_heat_rate_eqn(x, h(i), params), 3000);
    V4(i) = g_load_eqn(40.0, h(i), params);
    V5(i) = g_load_eqn(15.0, h(i), params);
end

min_data = min(V3, V4);
x2 = [min_data, fliplr(V1)]/1000;
inBetween = [h, fliplr(h)]/1000;

figure;
fill(x2, inBetween, 'y');
hold on;
plot(y(4,:)/1000, (y(1,:) - params.Rp)/1000, 'LineWidth', 1.5);
plot(V1/1000, h/1000);
plot(V2/1000, h/1000);
plot(V3/1000, h/1000, '--');
plot(V4/1000, h/1000, '--');
plot(V5/1000, h/1000, '--');
hold off;
xlim([0, V0/1000*1.025]);
xlabel('V [km/s]');
ylabel('h [km]');
name = legend('', 'Numerical Solution', 'Glide, $\alpha = 40^{\circ}$', ...
        'Glide, $\alpha = 10^{\circ}$', ...
        "$\dot{Q}_{max} = " + params.Qdot_max + "\;kW/m^2$", ...
        "$n_{max} = " + params.n_max + "$, $\alpha = 40^{\circ}$", ...
        "$n_{max} = " + params.n_max + "$, $\alpha = 10^{\circ}$");
set(name,'Interpreter','latex');
grid;

dfgdfg

%% Problem definition
problem.dynamics = @(t, x, u) dynamics(t, x, u, params);
problem.pathObj = @(t, x, u) pathObj(t, x, u, params);
problem.bndObj = @bndObj;
problem.pathCst = @(t, x, u) pathCst(t, x, u, params);

problem.bounds.t0.low = 0;
problem.bounds.t0.upp = 0;
problem.bounds.tf.low = 0;
problem.bounds.tf.upp = 600;

problem.bounds.x0.low = params.x0;
problem.bounds.x0.upp = params.x0;
problem.bounds.xf.low = [Rf; tauf; deltaf; Vf/2; -deg2rad(5.0); -pi/2];
problem.bounds.xf.upp = [Rf; tauf; deltaf; Vf*3; 0; pi/2];

problem.bounds.x.low = [Rf; tau0; delta0; Vf; -deg2rad(15.0); -pi/2];
problem.bounds.x.upp = [R0; tauf*1.5; deltaf*1.5; V0*1.5; deg2rad(0.0); pi/2];

problem.bounds.u.low = [deg2rad(-180)];
problem.bounds.u.upp = [deg2rad(180)];

problem.guess.t = ts;
problem.guess.x = y;
problem.guess.u = u;

problem.options(1).nIntervals   = 25;
problem.options(1).scheme       = 'trapezoid';
problem.options(1).solver       = 'fmincon';

problem.options(1).opt.Algorithm               = 'sqp';
problem.options(1).opt.EnableFeasibilityMode   = true;
problem.options(1).opt.SubproblemAlgorithm     = 'cg';
problem.options(1).opt.Display                 = 'iter';
problem.options(1).opt.MaxFunctionEvaluations  = 1e8;
problem.options(1).opt.MaxIterations           = 10000;
problem.options(1).opt.UseParallel             = true;

% problem.options(2).nIntervals   = 500;
% problem.options(2).scheme       = 'trapezoid';
% problem.options(2).solver       = 'fmincon';
% 
% problem.options(2).opt.Algorithm               = 'interior-point';
% problem.options(2).opt.Display                 = 'iter';
% problem.options(2).opt.MaxFunctionEvaluations  = 1e8;
% problem.options(2).opt.MaxIterations           = 10000;
% problem.options(2).opt.UseParallel             = true;

sol = direct_collocation(problem);

params.n_max = 4;
problem.pathCst = @(t, x, u) pathCst(t, x, u, params);
problem.guess.t = linspace(sol.t(1), sol.t(end), 500);
problem.guess.x = sol.interp.x(problem.guess.t);
problem.guess.u = sol.interp.u(problem.guess.t);
problem.options(1).opt.Algorithm               = 'sqp';
problem.options(1).nIntervals   = 25;
sol = direct_collocation(problem);

params.n_max = 3.5;
problem.pathCst = @(t, x, u) pathCst(t, x, u, params);
problem.guess.t = linspace(sol.t(1), sol.t(end), 500);
problem.guess.x = sol.interp.x(problem.guess.t);
problem.guess.u = sol.interp.u(problem.guess.t);
problem.options(1).opt.Algorithm               = 'sqp';
problem.options(1).nIntervals   = 25;
sol = direct_collocation(problem);

params.n_max = 3;
problem.pathCst = @(t, x, u) pathCst(t, x, u, params);
problem.guess.t = linspace(sol.t(1), sol.t(end), 500);
problem.guess.x = sol.interp.x(problem.guess.t);
problem.guess.u = sol.interp.u(problem.guess.t);
problem.options(1).opt.Algorithm               = 'sqp';
problem.options(1).nIntervals   = 25 ;
sol = direct_collocation(problem);

%% Plot results
N = 1000;
t = linspace(sol.t(1), sol.t(end), N);
x = sol.interp.x(t);
u = sol.interp.u(t);
err = sol.interp.err(t);
constrs = calc_constraints(t, x, u, params);

h = linspace(R0-params.Rp, Rf-params.Rp, 100);
V1 = zeros(1, 100);
V2 = zeros(1, 100);
V3 = zeros(1, 100);
V4 = zeros(1, 100);
V5 = zeros(1, 100);

for i=1:100
    V1(i) = fzero(@(x) glide_eqn(40.0, x, h(i), params), 3000);
    V2(i) = fzero(@(x) glide_eqn(15.0, x, h(i), params), 3000);
    V3(i) = fzero(@(x) total_heat_rate_eqn(x, h(i), params), 3000);
    V4(i) = g_load_eqn(40.0, h(i), params);
    V5(i) = g_load_eqn(15.0, h(i), params);
end

figure;
title('Position Residuals');
subplot(3, 1, 1);
plot(t, err(1,:)/1000);
xlabel('t [s]');
ylabel('R [km/s]');
subplot(3, 1, 2);
plot(t, rad2deg(err(2,:)));
xlabel('t [s]');
ylabel('\tau [deg/s]');
subplot(3, 1, 3);
plot(t, rad2deg(err(3,:)));
xlabel('t [s]');
ylabel('\delta [deg/s]');

figure;
title('Velocity Residuals');
subplot(3, 1, 1);
plot(t, err(4,:)/1000);
xlabel('t [s]');
ylabel('V [km/s^2]');
subplot(3, 1, 2);
plot(t, rad2deg(err(5,:)));
xlabel('t [s]');
ylabel('\gamma [deg/s]');
subplot(3, 1, 3);
plot(t, rad2deg(err(6,:)));
xlabel('t [s]');
ylabel('\chi [deg/s]');

figure;
subplot(2, 1, 1);
plot(t, (x(1,:) - params.Rp)/1000);
xlabel('t [s]');
ylabel('h [km]');
subplot(2, 1, 2);
plot(t, x(4,:)/1000);
xlabel('t [s]');
ylabel('V [km/s]');

figure;
subplot(2, 1, 1);
plot(t, rad2deg(x(2,:)));
xlabel('t [s]');
ylabel('lon [deg]');
subplot(2, 1, 2);
plot(t, rad2deg(x(3,:)));
xlabel('t [s]');
ylabel('lat [deg]');

figure;
subplot(2, 1, 1);
plot(t, rad2deg(x(5,:)));
xlabel('t [s]');
ylabel('\gamma [deg]');
subplot(2, 1, 2);
plot(t, rad2deg(x(6,:)));
xlabel('t [s]');
ylabel('\chi [deg]');

figure;
title('Control input u(t)');
plot(t, rad2deg(u(1,:)));
xlabel('t [s]');
ylabel('\sigma [deg]');

figure;
title('Constraints');
subplot(3, 1, 1);
plot(t, constrs(1,:));
xlabel('t [s]');
ylabel('q [kW/m^2]');
subplot(3, 1, 2);
plot(t, constrs(2,:));
xlabel('t [s]');
ylabel('q_{\infty} [Pa]');
subplot(3, 1, 3);
plot(t, constrs(3,:));
xlabel('t [s]');
ylabel('n_{load} [-]');

min_data = min(V3, V4);
x2 = [min_data, fliplr(V1)]/1000;
inBetween = [h, fliplr(h)]/1000;

figure;
fill(x2, inBetween, 'y');
hold on;
plot(x(4,:)/1000, (x(1,:) - params.Rp)/1000, 'LineWidth', 1.5);
plot(V1/1000, h/1000);
plot(V2/1000, h/1000);
plot(V3/1000, h/1000, '--');
plot(V4/1000, h/1000, '--');
plot(V5/1000, h/1000, '--');
hold off;
xlim([0, V0/1000*1.025]);
xlabel('V [km/s]');
ylabel('h [km]');
name = legend('', 'Numerical Solution', 'Glide, $\alpha = 40^{\circ}$', ...
        'Glide, $\alpha = 10^{\circ}$', ...
        "$\dot{Q}_{max} = " + params.Qdot_max + "\;kW/m^2$", ...
        "$n_{max} = " + params.n_max + "$, $\alpha = 40^{\circ}$", ...
        "$n_{max} = " + params.n_max + "$, $\alpha = 10^{\circ}$");
set(name,'Interpreter','latex');
grid;

%% Functions
function Jp = pathObj(t, x, u, p)
    Jp = zeros(size(t));
end

function Jb = bndObj(x0, t0, xf, tf)
    Jb = 0;%xf(4);
end

function [Ca, Cb] = pathCst(t, x, u, p)
    g = gravity(x(1,:));
    h = x(1,:) - p.Rp;
    T = p.T_interp(h);
    rho = p.rho_interp(h);
    M = x(4,:) ./ sqrt(p.gamma.*p.Rg.*T);
    AoA = AoA_curve(M);

    q_inf = 0.5.*rho.*x(4,:).^2;
    CD = p.drag_clean(AoA, M);
    CL = p.lift_clean(AoA, M);
    D = q_inf.*CD.*p.Sref;
    L = q_inf.*CL.*p.Sref;
    
    Qdot = 1.9027e-4 .* sqrt(rho./p.Rn) .* x(4,:).^3;   % [W/m^2]
    AoA_rad = deg2rad(AoA);
    n = (D.*sin(AoA_rad) + L.*cos(AoA_rad)) ./ (p.m.*p.g0);               % [-]

    Ca = [Qdot./1000 - p.Qdot_max; q_inf - p.q_max; n - p.n_max];
    Cb = [];
end

function C = calc_constraints(t, x, u, p)
    h = x(1,:) - p.Rp;
    T = p.T_interp(h);
    rho = p.rho_interp(h);
    M = x(4,:) ./ sqrt(p.gamma.*p.Rg.*T);
    AoA = AoA_curve(M);

    q_inf = 0.5.*rho.*x(4,:).^2;
    CD = p.drag_clean(AoA, M);
    CL = p.lift_clean(AoA, M);
    D = q_inf.*CD.*p.Sref;
    L = q_inf.*CL.*p.Sref;
    
    Qdot = 1.9027e-4 .* sqrt(rho./p.Rn) .* x(4,:).^3;     % [W/m^2]
    AoA_rad = deg2rad(AoA);
    n = (D.*sin(AoA_rad) + L.*cos(AoA_rad)) ./ (p.m.*p.g0);               % [-]

    C = [Qdot./1000; q_inf; n];
end

function [value, isterminal, direction] = myEvent(t, x, p)
    value      = x(1) - p.xf(1);
    isterminal = 1;   % Stop the integration
    direction  = -1;
end