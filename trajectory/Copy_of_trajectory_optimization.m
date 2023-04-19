clear
close all
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
params.q_max = 5000;                % [Pa]
params.n_max = 2.5;                 % [-]

params.drag_clean = drag_clean_interpolant('spline');
params.lift_clean = lift_clean_interpolant('spline');
[params.T_interp, params.p_interp, params.rho_interp] = atmosphere('spline');

%% Initial conditions
R0 = params.Rp + 120.0e3;           % [m]
tau0 = deg2rad(10.3);               % [rad]
delta0 = deg2rad(0);             % [rad]
V0 = 7000;                          % [m/s]
gamma0 = deg2rad(-10.5);             % [rad]
chi0 = deg2rad(33);               % [rad]
params.x0 = [R0, tau0, delta0, V0, gamma0, chi0]';

%% Terminal conditions
Rf = params.Rp + 10e3;               % [m]
tauf = deg2rad(21);               % [rad]
deltaf = deg2rad(20);             % [rad]
Vf = 800.0;                         % [m/s]
gammaf = deg2rad(-2.0);             % [rad]
chif = deg2rad(90.0);               % [rad]
params.xf = [Rf, tauf, deltaf, Vf, gammaf, chif]';

%% Calculate the initial guess
N = 1000;
tspan = [0, 500];
ts = linspace(tspan(1), tspan(2), N);
options = odeset('MaxStep',2, 'Events', @(t, x) myEvent(t, x, params));
u = zeros(1,N);
y = zeros(6,N);
y(:,1) = params.x0;

for i=2:N
    h = y(1,i-1) - params.Rp;

    if h <= 60e3
        T = params.T_interp(h);
        rho = params.rho_interp(h);
        g = gravity(y(1,i-1));
    
        M = y(4,i-1) ./ sqrt(params.gamma.*params.Rg.*T);
        AoA = AoA_curve(M);
    
        q_inf = 0.5.*rho.*y(4,i-1).^2;
        CL = params.lift_clean(AoA, M);
        L = q_inf.*CL.*params.Sref;
    
        var = (g - y(4,i-1).^2./y(1,i-1)) .* cos(y(5,i-1)) .* params.m ./ L;
        tmp = min(max(var, -1.0), 1.0);
        
        d_chi = delta_heading(y(:,i-1), params);
        u(1,i) = acos(tmp) .* sgn(d_chi, deg2rad(5), u(1,i-1));
    end

    % Solve the problem
    solution = ode45(@(t, x) dynamics_fixed(t, x, u(1,i), params), [ts(i-1), ts(i)], y(:,i-1), options);

    % Save solution
    y(:,i) = solution.y(:,end);

    if y(1,i) <= params.xf(1)
        ts = ts(1:i);
        y = y(:,1:i);
        u = u(1:i);
        break;
    end
end

figure();
plot(ts, rad2deg(u));
title('Control Law: u(t)');
xlabel('t [s]');
ylabel('Bank Angle [deg]');
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
legend(['chi', 'chi_d', 'delta_chi']);
grid;

figure;
plot(ts, rad2deg(delta_heading(y, params)));
title('Heading angle error delta_chi(t)');
xlabel('t [s]');
ylabel('Heading Angle Error [deg]');
grid;

%% Reentry zone
% M = 5;
% N = 50;
% V_range = linspace(5000, 7000, M);
% gamma_range = linspace(deg2rad(-20.0), deg2rad(-5.0), N);
% feasible = zeros(M, N);
% 
% for i=1:M
%     for j=1:N
%         fprintf('Problem i=%i, j=%i\n', i, j);
% 
%         %% Solve each problem using given initial conditions
%         params.x0 = [R0, tau0, delta0, V_range(i), gamma_range(j), chi0]';
% 
%         grid_sz = 200;
%         tspan = [0, 500];
%         ts = linspace(tspan(1), tspan(2), grid_sz);
%         options = odeset('MaxStep',2, 'Events', @(t, x) myEvent(t, x, params));
%         u = zeros(1,grid_sz);
%         y = zeros(6,grid_sz);
%         y(:,1) = params.x0;
%         
%         for k=2:grid_sz
%             h = y(1,k-1) - params.Rp;
%         
%             if h <= 60e3
%                 T = params.T_interp(h);
%                 rho = params.rho_interp(h);
%                 g = gravity(y(1,k-1));
%             
%                 M = y(4,k-1) ./ sqrt(params.gamma.*params.Rg.*T);
%                 AoA = AoA_curve(M);
%             
%                 q_inf = 0.5.*rho.*y(4,k-1).^2;
%                 CL = params.lift_clean(AoA, M);
%                 L = q_inf.*CL.*params.Sref;
%             
%                 var = (g - y(4,k-1).^2./y(1,k-1)) .* cos(y(5,k-1)) .* params.m ./ L;
%                 tmp = min(max(var, -1.0), 1.0);
%                 
%                 d_chi = delta_heading(y(:,k-1), params);
%                 u(1,k) = acos(tmp) .* sgn(d_chi, deg2rad(5), u(1,k-1));
%             end
%         
%             % Solve the problem
%             solution = ode45(@(t, x) dynamics_fixed(t, x, u(1,k), params), [ts(k-1), ts(k)], y(:,k-1), options);
%         
%             % Save solution
%             y(:,k) = solution.y(:,end);
%         
%             if y(1,k) <= params.xf(1)
%                 ts = ts(1:k);
%                 y = y(:,1:k);
%                 u = u(1:k);
% 
%                 % Mark V-gamma pair as feasible
%                 feasible(i,j) = 1;
%                 break;
%             end
%         end
%     end
% end
% 
% figure;
% colorMap = [1 1 1; 0 1 0];
% colormap(colorMap);
% imagesc(rad2deg(gamma_range), V_range, feasible);

%gfkf
%% Problem definition
problem.func.dynamics = @(t, x, u) dynamics(t, x, u, params);
problem.func.pathObj = @(t, x, u) pathObj(t, x, u, params);
problem.func.bndObj = @bndObj;
problem.func.pathCst = @(t, x, u) pathCst(t, x, u, params);

problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = 0;
problem.bounds.finalTime.upp = Inf;

problem.bounds.initialState.low = params.x0;
problem.bounds.initialState.upp = params.x0;

problem.bounds.state.low = [Rf; tau0; delta0; Vf; -deg2rad(15.0); -pi/2];
problem.bounds.state.upp = [R0; tauf*1.5; deltaf*1.5; V0*1.5; deg2rad(0.0); pi/2];

problem.bounds.finalState.low = [Rf; tauf; deltaf; Vf/2; -deg2rad(5.0); -pi/2];
problem.bounds.finalState.upp = [Rf; tauf; deltaf; Vf*3; 0; pi/2];

problem.bounds.control.low = [-deg2rad(0.0)];
problem.bounds.control.upp = [deg2rad(180)];

problem.guess.time = ts;
problem.guess.state = y;
problem.guess.control = abs(u);
% problem.guess.time = [0, 55, 300];
% problem.guess.state = [R0, params.Rp + 7e4, Rf;
%                        tau0, (tau0+tauf)/2, tauf;
%                        delta0, (delta0+deltaf)/2, deltaf;
%                        V0, V0, Vf*2;
%                        gamma0, gamma0, gammaf;
%                        chi0, chi0, chi0];
% problem.guess.control = [deg2rad(65), deg2rad(65), deg2rad(65)];

% Solver procedure selection
% problem.options(1).verbose = 3;
% problem.options(1).method = 'rungeKutta';
% problem.options(1).rungeKutta.nSegment = 10;
% problem.options(1).rungeKutta.nSubStep = 2;
% 
% problem.options(1).verbose = 3;
% problem.options(1).method = 'rungeKutta';
% problem.options(1).rungeKutta.nSegment = 30;
% problem.options(1).rungeKutta.nSubStep = 1;
% problem.options(1).nlpOpt.Algorithm = 'sqp';
% problem.options(1).nlpOpt.MaxIter = 50;

problem.options(1).verbose = 3;
problem.options(1).method = 'rungeKutta';
problem.options(1).rungeKutta.nSegment = 100;
problem.options(1).rungeKutta.nSubStep = 1;
problem.options(1).nlpOpt.MaxIter = 50;

% problem.options(3).verbose = 3;
% problem.options(3).method = 'rungeKutta';
% problem.options(3).rungeKutta.nSegment = 100;
% problem.options(3).rungeKutta.nSubStep = 1;
% problem.options(3).nlpOpt.Algorithm = 'sqp';

%problem.options(1).nlpOpt.MaxIter = 1;
%problem.options(1).nlpOpt.Algorithm = 'interior-point';
%problem.options(1).nlpOpt.EnableFeasibilityMode = true;
%problem.options(1).nlpOpt.SubproblemAlgorithm = 'cg';

% problem.options(2).verbose = 3;
% problem.options(2).method = 'rungeKutta';
% problem.options(2).rungeKutta.nSegment = 100;
% problem.options(2).rungeKutta.nSubStep = 1;

% problem.options(1).verbose = 3;
% problem.options(1).method = 'rungeKutta';
% problem.options(1).rungeKutta.nSegment = 100;
% problem.options(1).rungeKutta.nSubStep = 1;

% problem.options(4).verbose = 3;
% problem.options(4).method = 'rungeKutta';
% problem.options(4).rungeKutta.nSegment = 200;
% problem.options(4).rungeKutta.nSubStep = 1;

% problem.options(5).verbose = 3;
% problem.options(5).method = 'rungeKutta';
% problem.options(5).rungeKutta.nSegment = 500;
% problem.options(5).rungeKutta.nSubStep = 1;

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
x = soln(end).interp.state(t);
u = soln(end).interp.control(t);
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
    V1(i) = fzero(@(x) glide_eqn(45.0, x, h(i), params), 3000);
    V2(i) = fzero(@(x) glide_eqn(15.0, x, h(i), params), 3000);
    V3(i) = fzero(@(x) total_heat_rate_eqn(x, h(i), params), 3000);
    V4(i) = g_load_eqn(45.0, h(i), params);
    V5(i) = g_load_eqn(15.0, h(i), params);
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
plot(t, rad2deg(u(1,:)) .* sgn(delta_heading(x, params), deg2rad(5), u));
xlabel('t [s]');
ylabel('\sigma [deg]');

% figure;
% title('Heading angle error delta_chi(t)');
% plot(t, rad2deg(delta_chi(x, params)));
% xlabel('t [s]');
% ylabel('\delta_chi [deg]');
% 
% figure;
% title('banking angle sign, sgn(u(t))');
% plot(t, sgn(delta_chi(x, params), deg2rad(20)));
% xlabel('t [s]');
% ylabel('sgn(\sigma) [-]');

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
%     g = gravity(x(1,:));
% 
%     h = x(1,:) - p.Rp;
%     T = p.T_interp(h);
%     rho = p.rho_interp(h);
%     M = x(4,:) ./ sqrt(p.gamma.*p.Rg.*T);
%     AoA = AoA_curve(M);
% 
%     d_chi = delta_heading(x, p);
%     sigma = u(1,:) .* sgn(d_chi, deg2rad(5), u(1,:));
% 
%     V_squared = x(4,:).^2;
%     q_inf = 0.5.*rho.*V_squared;
%     CL = p.lift_clean(AoA, M);
%     CD = p.drag_clean(AoA, M);
%     L = q_inf.*CL.*p.Sref;
%     D = q_inf.*CD.*p.Sref;
%     
%     Jp = -D + L .* cos(u(1,:));
% 
%     Jp =  (L./p.m.*cos(sigma) + (V_squared./x(1,:) - g) .* cos(x(5,:))).^2;
    Jp = zeros(size(t));
end

function Jb = bndObj(t0, x0, tf, xf)
    Jb = xf(4);
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
    
    Qdot = 1.9027e-4 .* sqrt(rho./p.Rn) .* x(4,:).^3;     % [W/m^2]
    n = sqrt(D.^2 + L.^2) ./ (p.m.*p.g0);               % [-]
    %Vg_dot = L./p.m.*cos(u(1,:)) + (x(4,:).^2./x(1,:) - g) .* cos(x(5,:));

    %AoA = ones(size(M)) * 10;
    %CL_g = p.lift_clean(AoA, M);
    %CD_g = p.drag_clean(AoA, M);

    %g_constr = rho.*x(4,:).^2.*p.Sref ./ (2.*p.m.*p.g0) .* sqrt(CL_g.^2 + CD_g.^2) - p.n_max;
    u_constr = u(1,:) .* (h > 6e4);

    Ca = [Qdot./1000 - p.Qdot_max; q_inf - p.q_max; n - p.n_max];
    %Ca = [];
    Cb = [u_constr];
end

function C = calc_constraints(t, x, u, p)
    h = x(1,:) - p.Rp;
    T = p.T_interp(h);
    rho = p.rho_interp(h);
    M = x(4,:) ./ sqrt(p.gamma.*p.Rg.*T);
    AoA = ones(size(M)) * 40.0;

    q_inf = 0.5.*rho.*x(4,:).^2;
    CD = p.drag_clean(AoA, M);
    CL = p.lift_clean(AoA, M);
    D = q_inf.*CD.*p.Sref;
    L = q_inf.*CL.*p.Sref;
    
    Qdot = 1.9027e-4 .* sqrt(rho./p.Rn) .* x(4,:).^3;     % [W/m^2]
    n = sqrt(D.^2 + L.^2) ./ (p.m.*p.g0);               % [-]

    C = [Qdot./1000; q_inf; n];
end

function [value, isterminal, direction] = myEvent(t, x, p)
    value      = x(1) - p.xf(1);
    isterminal = 1;   % Stop the integration
    direction  = -1;
end