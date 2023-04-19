clear
close all
addpath(pwd + "\fcn\aerodynamic_coefficients")
addpath(pwd + "\fcn\atmospheric_model")
addpath(pwd + "\fcn\gravity_model")

%% Parameters setup
params.mu = 0.042828e15;            % [m^3/s^2]
params.Rp = 3396.2 * 1e3;           % [m]
%params.omega_cb = 7.07765809e-5;    % [rad/s]
params.omega_cb = 0;
params.gamma = 1.2941;              % [-]
params.Rg = 8314.3/44.01;           % [-]
params.Sref = 110;                  % [m^2]
params.Rn = 1.2;                    % [m]
params.m = 26029.0;                 % [kg]
params.g0 = 9.81;                   % [m/s^2]

params.Qdot_max = 50;               % [kW/m^2]
params.q_max = 5000;                % [Pa]
params.n_max = 1.5;                 % [-]

%% Initial conditions
R0 = params.Rp + 120.0e3;           % [m]
tau0 = deg2rad(10.3);               % [rad]
delta0 = deg2rad(0);             % [rad]
%V0 = sqrt(params.mu / R0);         % [m/s]
V0 = 3800;                          % [m/s]
gamma0 = deg2rad(-1.0);             % [rad]
chi0 = deg2rad(33);               % [rad]

%% Terminal conditions
Rf = params.Rp + 10e3;               % [m]
tauf = deg2rad(90.3);               % [rad]
deltaf = deg2rad(56);             % [rad]
Vf = 700.0;                         % [m/s]
gammaf = deg2rad(-1.0);             % [rad]
chif = deg2rad(90.0);               % [rad]

%% Characteristic variables
R_s = R0;
tau_s = pi;
delta_s = deg2rad(89);
V_s = V0;
gamma_s = deg2rad(10);
chi_s = pi;
sigma_s = deg2rad(90);

params.R_s = R_s;
params.tau_s = tau_s;
params.delta_s = delta_s;
params.V_s = V_s;
params.gamma_s = gamma_s;
params.chi_s = chi_s;
params.sigma_s = sigma_s;

%% Problem definition
problem.func.dynamics = @(t, x, u) dynamics(t, x, u, params);
problem.func.pathObj = @(t, x, u) pathObj(t, x, u, params);
problem.func.bndObj = @bndObj;
problem.func.pathCst = @(t, x, u) pathCst(t, x, u, params);

problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = 0;
problem.bounds.finalTime.upp = Inf;

problem.bounds.initialState.low = [R0/R_s; tau0/tau_s; delta0/delta_s; V0/V_s; gamma0/gamma_s; chi0/chi_s];
problem.bounds.initialState.upp = [R0/R_s; tau0/tau_s; delta0/delta_s; V0/V_s; gamma0/gamma_s; chi0/chi_s];

problem.bounds.state.low = [Rf/R_s; -1; -1;     Vf/V_s; -1; -1];
problem.bounds.state.upp = [R0/R_s;  1;  1; V0/V_s*1.5;  1;  1];

problem.bounds.finalState.low = [Rf/R_s; tauf/tau_s; deltaf/delta_s; Vf/V_s; -1; chif/chi_s];
problem.bounds.finalState.upp = [Rf/R_s; tauf/tau_s; deltaf/delta_s; Vf/V_s;  0; chif/chi_s];

problem.bounds.control.low = [-1.0];
problem.bounds.control.upp = [1.0];

problem.guess.time = [0, 5000];
problem.guess.state = [        R0/R_s,    Rf/R_s;
                           tau0/tau_s,       0.0;
                       delta0/delta_s,       0.0;
                               V0/V_s,    Vf/V_s;
                       gamma0/gamma_s,       0.0;
                           chi0/chi_s,       0.0];

problem.guess.control = [0.0, 0.0];

%problem.options.nlpOpt.options.MaxIterations = 1000;

% Solver procedure selection
problem.options(1).verbose = 3;
problem.options(1).method = 'rungeKutta';
problem.options(1).rungeKutta.nSegment = 10;
problem.options(1).rungeKutta.nSubStep = 2;
problem.options(1).nlpOpt.MaxIterations = 1000;
%problem.options(1).nlpOpt.EnableFeasibilityMode = true;
%problem.options(1).chebyshev.nColPts = 10;

problem.options(2).verbose = 3;
problem.options(2).method = 'hermiteSimpson';
problem.options(2).hermiteSimpson.nSegment = 50;
%problem.options(2).chebyshev.nColPts = 50;

% problem.options(3).verbose = 3;
% problem.options(3).method = 'chebyshev';
% problem.options(3).chebyshev.nColPts = 16;
% 
% problem.options(4).verbose = 3;
% problem.options(4).method = 'chebyshev';
% problem.options(4).chebyshev.nColPts = 20;
% 
% problem.options(5).verbose = 3;
% problem.options(5).method = 'chebyshev';
% problem.options(5).chebyshev.nColPts = 30;
% 
% problem.options(6).verbose = 3;
% problem.options(6).method = 'chebyshev';
% problem.options(6).chebyshev.nColPts = 100;
% problem.options(6).defaultAccuracy = 'high';

soln = optimTraj(problem);

%% Plot results
N = 1000;
t = linspace(0, soln(end).grid.time(end), N);
x = soln(end).interp.state(t) .* [R_s, tau_s, delta_s, V_s, gamma_s, chi_s]';
u = soln(end).interp.control(t) .* [sigma_s]';
constrs = calc_constraints(t, x, u, params);

dx_dt = diff(x,1,2) ./ diff(t);
real_dx_dt = dynamics(t, x, u, params);
diff_dx_dt = abs(dx_dt - real_dx_dt(:,1:end-1));

h = linspace(R0-params.Rp, Rf-params.Rp, 100);
V1 = zeros(1, 100);
V2 = zeros(1, 100);
V3 = zeros(1, 100);
V4 = zeros(1, 100);
V5 = zeros(1, 100);

for i=1:100
    V1(i) = fzero(@(x) glide_eqn(40.0, x, h(i), params), 3000);
    V2(i) = fzero(@(x) glide_eqn(10.0, x, h(i), params), 3000);
    V3(i) = fzero(@(x) total_heat_rate_eqn(x, h(i), params), 3000);
    V4(i) = g_load_eqn(40.0, h(i), params);
    V5(i) = g_load_eqn(10.0, h(i), params);
end

figure;
title('Position Residuals');
subplot(3, 1, 1);
plot(t(1:end-1), diff_dx_dt(1,:)/1000);
xlabel('t [s]');
ylabel('h [km/s]');
subplot(3, 1, 2);
plot(t(1:end-1), rad2deg(diff_dx_dt(2,:)));
xlabel('t [s]');
ylabel('\tau [deg/s]');
subplot(3, 1, 3);
plot(t(1:end-1), rad2deg(diff_dx_dt(3,:)));
xlabel('t [s]');
ylabel('\delta [deg/s]');

figure;
title('Velocity Residuals');
subplot(3, 1, 1);
plot(t(1:end-1), diff_dx_dt(4,:)/1000);
xlabel('t [s]');
ylabel('V [km/s^2]');
subplot(3, 1, 2);
plot(t(1:end-1), rad2deg(diff_dx_dt(5,:)));
xlabel('t [s]');
ylabel('\gamma [deg/s]');
subplot(3, 1, 3);
plot(t(1:end-1), rad2deg(diff_dx_dt(6,:)));
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
name = legend('', 'Solution', 'Glide, $\alpha = 40^{\circ}$', ...
        'Glide, $\alpha = 10^{\circ}$', ...
        "$\dot{Q}_{max} = " + params.Qdot_max + "\;kW/m^2$", ...
        "$n_{max} = " + params.n_max + "$, $\alpha = 40^{\circ}$", ...
        "$n_{max} = " + params.n_max + "$, $\alpha = 10^{\circ}$");
set(name,'Interpreter','latex');
grid;

%% Functions
function g = gravity(R)
    mu = 0.042828e15;    % [m^3/s^2]
    g = mu./R.^2;
end

function res = AoA_curve(M)
    f = 15;
    g = 40;
    
    K = 1;
    M_c = 7;

    s = (1 + tanh(K.*(M - M_c))) ./ 2;
    res = deg2rad(g.*s + (1-s).*f);
end

function res = glide_eqn(AoA, V, h, p)
    R = p.Rp + h;
    g = gravity(R);
    W = p.m .* g;
    [T, ~, rho] = mars_atmosphere(h);
    M = V ./ sqrt(p.gamma .* p.Rg .* T);
    CL = lift_clean(AoA, M);

    res = 2.*(W./p.Sref)./CL .* (1./V.^2 - 1./(R.*g)) - rho;
end

function res = g_load_eqn(AoA, h, p)
    [T, ~, rho] = mars_atmosphere(h);
    M = 2000 ./ sqrt(p.gamma .* p.Rg .* T);
    CL = lift_clean(AoA, M);
    CD = drag_clean(AoA, M);

    res = sqrt(2.*p.m.*p.g0.*p.n_max ./ (rho .* p.Sref * sqrt(CL.^2 + CD.^2)));
end

function res = total_heat_rate_eqn(V, h, p)
    [~, ~, rho] = mars_atmosphere(h);
    q_conv = 1.9027e-4 .* sqrt(rho./p.Rn) .* V.^3;
    res = q_conv ./ 1000 - p.Qdot_max;
end

% function [rho, T] = density(h)
%     if h < 7000
%         T = -31 - 0.000998.*h;
%     else
%         T = -23.4 - 0.00222.*h;
%     end
% 
%     p = 0.699 .* exp(-0.00009.*h);
%     rho = p ./ (0.1921 .* (T + 273.1));
% end

function xdot = dynamics(t, x, u, p)
    % x = [R, lon, lat, V, gamma, chi]
    % u = [AoA, sigma]
    
    sd = sin(x(3,:) .* p.delta_s);
    cd = cos(x(3,:) .* p.delta_s);
    td = tan(x(3,:) .* p.delta_s);
    sg = sin(x(5,:) .* p.gamma_s);
    cg = cos(x(5,:) .* p.gamma_s);
    sc = sin(x(6,:) .* p.chi_s);
    cc = cos(x(6,:) .* p.chi_s);
    
    R = x(1,:) .* p.R_s;
    V = x(4,:) .* p.V_s;
    sigma = u(1,:) .* p.sigma_s;

    h = R - p.Rp;
    [T, ~, rho] = mars_atmosphere(h);
    g = gravity(R);
    
    M = V ./ sqrt(p.gamma.*p.Rg.*T);
    AoA = 40.0; % Fixed

    q_inf = 0.5.*rho.*V.^2;
    CD = drag_clean(AoA, M);
    CL = lift_clean(AoA, M);
    D = q_inf.*CD.*p.Sref;
    L = q_inf.*CL.*p.Sref;

    Fv = -D - p.m.*g.*sg;
    Fg = L.*cos(sigma) - p.m.*g.*cg;
    Fc = -L.*sin(sigma);
    
    xdot = zeros(size(x));
    xdot(1,:) = V.*sg ./ p.R_s;
    xdot(2,:) = V.*sc.*cg ./ (R.*cd) ./ p.tau_s;
    xdot(3,:) = V.*cc.*cg ./ R ./ p.delta_s;
    xdot(4,:) = (Fv./p.m + p.omega_cb.^2.*R.*cd.*(sg.*cd - cg.*sd.*cc)) ./ p.V_s;
    xdot(5,:) = (Fg./p.m + 2.*p.omega_cb.*V.*cd.*sc + V.^2./R.*cg + ...
        p.omega_cb.^2.*R.*cd.*(cd.*cg + sg.*sd.*cc))./V ./ p.gamma_s;
    xdot(6,:) = (Fc./p.m + 2.*p.omega_cb.*V.*(sd.*cg - cd.*sg.*cc) + ...
        V.^2./R.*cg.^2.*td.*sc + p.omega_cb.^2.*R.*cd.*sd.*sc)./(V.*cg) ./ p.chi_s;
end

function Jp = pathObj(t, x, u, p)
    R = x(1,:) .* p.R_s;
    V = x(4,:) .* p.V_s;
    gamma = x(5,:) .* p.gamma_s;
    sigma = u(1,:) .* p.sigma_s;
    g = gravity(R);

    h = R - p.Rp;
    [T, ~, rho] = mars_atmosphere(h);
    M = V ./ sqrt(p.gamma.*p.Rg.*T);
    AoA = 40.0;

    q_inf = 0.5.*rho.*V.^2;
    CL = lift_clean(AoA, M);
    L = q_inf.*CL.*p.Sref;

    Jp =  (L./p.m.*cos(sigma) + (V.^2./R - g) .* cos(gamma)).^2;
end

function Jb = bndObj(t0, x0, tf, xf)
    Jb = 0;
end

function [Ca, Cb] = pathCst(t, x, u, p)
    R = x(1,:) .* p.R_s;
    V = x(4,:) .* p.V_s;
    gamma = x(5,:) .* p.gamma_s;
    sigma = u(1,:) .* p.sigma_s;
    g = gravity(R);

    h = R - p.Rp;
    [T, ~, rho] = mars_atmosphere(h);
    M = V ./ sqrt(p.gamma.*p.Rg.*T);
    AoA = 40.0;

    q_inf = 0.5.*rho.*V.^2;
    CD = drag_clean(AoA, M);
    CL = lift_clean(AoA, M);
    D = q_inf.*CD.*p.Sref;
    L = q_inf.*CL.*p.Sref;
    
    Qdot = 1.9027e-4 .* sqrt(rho./p.Rn) .* V.^3;     % [W/m^2]
    n = sqrt(D.^2 + L.^2) ./ (p.m.*p.g0);               % [-]
    
    % Constraints on banking angle
    %sigma_dot = diff(rad2deg(sigma)) ./ diff(t);
    %sigma_dot = [sigma_dot, 0];
    %sigma_dot_max = 3.0;

    %CL_g = lift_clean(ones(size(AoA)) .* 10.0, M);
    %CD_g = drag_clean(ones(size(AoA)) .* 10.0, M);

    %g_constr = rho.*V.^2.*p.Sref ./ (2.*p.m.*p.g0) .* sqrt(CL_g.^2 + CD_g.^2) - p.n_max;
    gamma_dot = L./p.m.*cos(sigma) + (V.^2./R - g) .* cos(gamma);
    
    Ca = [Qdot./1000 - p.Qdot_max; q_inf - p.q_max; n - p.n_max; gamma_dot];
    Cb = [];
end

function C = calc_constraints(t, x, u, p)
    R = x(1,:);
    V = x(4,:);
    sigma = u(1,:);

    h = R - p.Rp;
    [T, ~, rho] = mars_atmosphere(h);
    M = V ./ sqrt(p.gamma.*p.Rg.*T);
    AoA = 40.0;

    q_inf = 0.5.*rho.*V.^2;
    CD = drag_clean(AoA, M);
    CL = lift_clean(AoA, M);
    D = q_inf.*CD.*p.Sref;
    L = q_inf.*CL.*p.Sref;
    
    Qdot = 1.9027e-4 .* sqrt(rho./p.Rn) .* V.^3;     % [W/m^2]
    n = sqrt(D.^2 + L.^2) ./ (p.m.*p.g0);               % [-]
    
    % Constraints on angle of attack and banking angle
    sigma_dot = diff(rad2deg(sigma)) ./ diff(t);
    sigma_dot = [sigma_dot, 0];

    C = [Qdot./1000; q_inf; n; sigma_dot];
end