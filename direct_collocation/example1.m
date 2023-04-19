clear;
close all;

problem.bndObj              = @bndObj;
problem.pathObj             = @pathObj;
problem.pathCst             = @pathCst;
problem.dynamics            = @dynamics;

problem.bounds.t0.low       = 0;
problem.bounds.t0.upp       = 0;
problem.bounds.tf.low       = 1;
problem.bounds.tf.upp       = 1;

problem.bounds.x0.low       = [0; 0];
problem.bounds.x0.upp       = [0; 0];
problem.bounds.xf.low       = [1; 0];
problem.bounds.xf.upp       = [1; 0];

problem.bounds.x.low        = [-1; -100];
problem.bounds.x.upp        = [1; 100];
problem.bounds.u.low        = -10;
problem.bounds.u.upp        = 10;

problem.guess.t             = [0, 1];
problem.guess.x             = [0, 1; 1 1];
problem.guess.u             = [0, 0];

opt.Algorithm = 'interior-point';
opt.Display = 'iter';
opt.MaxFunctionEvaluations = 1e6;

problem.options(1).nIntervals   = 10;
problem.options(1).scheme       = 'trapezoid';
problem.options(1).solver       = 'fmincon';
problem.options(1).opt          = opt;

problem.options(2).nIntervals   = 100;
problem.options(2).scheme       = 'simpson';
problem.options(2).solver       = 'fmincon';
problem.options(2).opt          = opt;

sol = direct_collocation(problem);

exact_sol_x = @(t) 3.*t.^2 - 2.*t.^3;
exact_sol_u = @(t) 6 - 12.*t;

t = linspace(sol.t(1), sol.t(end), 10000);
x = sol.interp.x(t);
u = sol.interp.u(t);
err = sol.interp.err(t);

figure;
plot(t, x(1,:));
hold on;
plot(t, exact_sol_x(t));
hold off;
xlabel('Time [s]');
ylabel('Position [m]');

figure;
plot(t, u);
hold on;
plot(t, exact_sol_u(t));
hold off;
xlabel('Time [s]');
ylabel('Force (u)');

figure;
plot(t, err(1,:));
xlabel('Time [s]');
ylabel('Position dynamics error [m/s]');

%% Problem-Specific Functions
function f = bndObj(x0, t0, xf, tf)
    f = 0;
end

function f = pathObj(t, x, u)
    f = u(1,:).^2;
end

function y = dynamics(t, x, u)
    y = [x(2,:); u(1,:)];
end

function [Cin, Ceq] = pathCst(t, x, u)
    Cin = [];
    Ceq = [];
end