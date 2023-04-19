x0 = simdb.nav.x0;
P0 = simdb.nav.P0;
Q = simdb.nav.Q;
R = simdb.nav.R;
dt = simdb.nav.tsamp
dq_K = simdb.nav.dq_K;
mu = simdb.env.mu;
Ma = simdb.nav.Ma;
Mo = simdb.nav.Mo;

tspan = [0, 100];
t = tspan(1);

x = x0;
P = P0;
idx = 0;

while t < tspan(2)
    % Predict
    x_est = fcn_est_vec(x, dt, dq_K, mu);
    F = fcn_est_matrix(x, dt, dq_K, mu);
    P_est = F*P + P*F' + Q;

    z_meas = zeros(6,1);
    z_est = fcn_meas_vec(x_est, dt, Ma, Mo);
    H = fcn_meas_matrix(x_est, dt, Ma, Mo);
    
    % Update
    K = P_est*H'*inv(H*P_est*H' + R);
    x = x + K*(z_meas - z_est);
    P = (eye(22) - K*H)*P_est;

    t = t + dt;
%     idx = idx + 1;
%     if idx >= 10
%         t
%         idx = 0;
%         x
%     end
end

x
P

function x_est = fcn_est_vec(x, dt, K, mu)
    x_est = zeros(size(x));

    g = -mu/norm(x(1:3))^3 * x(1:3);
    eps = 1 - (x(13)^2 + x(14)^2 + x(15)^2 + x(16)^2);

    x_est(1:3) = x(1:3) + x(4:6)*dt + 0.5*(x(7:9) + g)*dt^2;
    x_est(4:6) = x(4:6) + (x(7:9) + g)*dt;
    x_est(7:9) = x(7:9);
    x_est(10:12) = x(10:12);
    x_est(13:16) = x(13:16) + (0.5*quatmul([0; x(10:12)], x(13:16)) + K*eps*x(13:16))*dt;
    x_est(17:19) = x(17:19);
    x_est(20:22) = x(20:22);
end

function z = fcn_meas_vec(x, dt, Ma, Mo)
    z = zeros(6,1);
    I = eye(3);

    z(1:3) = ((I + Ma)*quatrot(quatconj(x(13:16)), x(7:9)) + x(17:19))*dt;
    z(4:6) = ((I + Mo)*x(10:12) + x(20:22))*dt;
end

function F = fcn_est_matrix(x, dt, K, mu)
    F = zeros(22);

    R = norm(x(1:3));
    I = eye(3);
    dg = mu/R^5*(3*x(1:3)*x(1:3)' - R^2*I);

    dq_dw = 0.5*[-x(14) -x(15) -x(16);
                  x(13)  x(16) -x(15);
                 -x(16)  x(13)  x(14);
                  x(15) -x(14)  x(13)];

    Omega = [   0  -x(10) -x(11) -x(12);
             x(10)     0  -x(12)  x(11);
             x(11)  x(12)     0  -x(10);
             x(12) -x(11)  x(10)     0];

    F(1:3,4:6) = I;

    F(4:6,1:3) = dg;
    F(4:6,7:9) = I;

    F(13:16,10:12) = dq_dw;
    F(13:16,13:16) = 0.5*Omega;
end

function H = fcn_meas_matrix(x, dt, Ma, Mo)
    H = zeros(6,22);
    
    I = eye(3);

    H(1:3,7:9) = (I + Ma) * [    1 - 2*(x(15)^2 + x(16)^2)       2*(x(14)*x(15) + x(13)*x(16))       2*(x(14)*x(16) - x(13)*x(15));
                             2*(x(14)*x(15) - x(13)*x(16))           1 - 2*(x(14)^2 + x(16)^2)       2*(x(15)*x(16) + x(13)*x(14));
                             2*(x(14)*x(16) + x(13)*x(15))       2*(x(15)*x(16) - x(13)*x(14))           1 - 2*(x(14)^2 + x(15)^2)] * dt;

    qw_x = x(13)*x(7);  qw_y = x(13)*x(8);  qw_z = x(13)*x(9);
    qx_x = x(14)*x(7);  qx_y = x(14)*x(8);  qx_z = x(14)*x(9);
    qy_x = x(15)*x(7);  qy_y = x(15)*x(8);  qy_z = x(15)*x(9);
    qz_x = x(16)*x(7);  qz_y = x(16)*x(8);  qz_z = x(16)*x(9);
    
    H(1:3,13:16) = 2*(I + Ma) * [ qz_y - qy_z             qy_y + qz_z    -2*qy_x + qx_y - qw_z    -2*qz_x + qw_y + qx_z;
                                 -qz_x + qx_z    qy_x - 2*qx_y + qw_z              qx_x + qz_z    -qw_x - 2*qz_y + qy_z;
                                  qy_x - qx_y    qz_x - qw_y - 2*qx_z     qw_x + qz_y - 2*qy_z              qx_x + qy_y] * dt;
    H(1:3,17:19) = I*dt;

    H(4:6,10:12) = (I + Mo)*dt;
    H(4:6,20:22) = I*dt;
end