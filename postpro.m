files = dir('output');
subFolders = files([files.isdir]);
outputFolder = subFolders(end).name;
workingDir = fullfile('output', outputFolder);

files = dir(workingDir);
subFolders = files([files.isdir]);
subFolderNames = {subFolders(3:end).name};
N = length(subFolderNames);

results = struct();
for i=1:N
    simFile = fullfile(workingDir, subFolderNames{i}, 'results.mat');
    inputFile = fullfile(workingDir, subFolderNames{i}, 'input.mat');
    results.(strcat("sim", string(i))) = load(simFile, 'out');
    results.(strcat("input", string(i))) = load(inputFile, 'simdb');
end

%% Instantiate figures
f_h      = figure(1);
f_dh     = figure(2);
f_V      = figure(3);
f_dV     = figure(4);
f_land   = figure(5);
f_rho    = figure(6);
f_alpha  = figure(7);
f_beta   = figure(8);
f_sigma  = figure(9);
f_nav_R  = figure(10);
f_nav_V  = figure(11);
f_nav_ae = figure(12);

%% Preprocessing
max_sz = 1;
for i=1:N
    sim = results.(['sim',num2str(i)]).out;
    simdb = results.(['input',num2str(i)]).simdb;
    ts = sim.sc_pos_vrt.Time;
    max_sz = max(max_sz, length(ts));
end

%% Create data arrays
taem        = zeros(N,2);
height_err  = zeros(max_sz,N);
vel_err     = zeros(max_sz,N);
dtau_err    = zeros(max_sz,N);
ddelta_err  = zeros(max_sz,N);
dgamma_err  = zeros(max_sz,N);
dchi_err    = zeros(max_sz,N);

nav_dR      = zeros(max_sz,N);
nav_dtau    = zeros(max_sz,N);
nav_ddelta  = zeros(max_sz,N);
nav_dV      = zeros(max_sz,N);
nav_dgamma  = zeros(max_sz,N);
nav_dchi    = zeros(max_sz,N);
nav_dalpha  = zeros(max_sz,N);
nav_dbeta   = zeros(max_sz,N);
nav_dsigma  = zeros(max_sz,N);
max_ts      = 0:0.05:(max_sz-1)/20;

%% Plot to figures
for i=1:N
    sim = results.(['sim',num2str(i)]).out;
    simdb = results.(['input',num2str(i)]).simdb;

    % Plot Height
    ts = sim.sc_pos_vrt.Time;
    data = sim.sc_pos_vrt.Data;
    set(0, "CurrentFigure", f_h);
    plot(ts, (data(:,1) - simdb.env.Rp)*1e-3, 'b');
    hold on;

    % Plot Height error
    data = sim.sc_pos_vrt.Data - sim.y_ref.Data(:,1:3);
    sz = length(data(:,1));
    height_err(1:sz,i) = data(:,1);
    dtau_err(1:sz,i) = rad2deg(data(:,2));
    ddelta_err(1:sz,i) = rad2deg(data(:,3));
    set(0, "CurrentFigure", f_dh);
    subplot(3,1,1);
    plot(ts, data(:,1), 'b');
    hold on;
    subplot(3,1,2);
    plot(ts, rad2deg(data(:,2)), 'b');
    hold on;
    subplot(3,1,3);
    plot(ts, rad2deg(data(:,3)), 'b');
    hold on;

    % Plot Velocity
    ts = sim.sc_vel_tg.Time;
    data = sim.sc_vel_tg.Data;
    set(0, "CurrentFigure", f_V);
    plot(ts, data(:,1), 'b');
    hold on;

    % Plot Velocity error
    data = sim.sc_vel_tg.Data - sim.y_ref.Data(:,4:6);
    sz = length(data(:,1));
    vel_err(1:sz,i) = data(:,1);
    dgamma_err(1:sz,i) = rad2deg(data(:,2));
    dchi_err(1:sz,i) = rad2deg(data(:,3));
    set(0, "CurrentFigure", f_dV);
    subplot(3,1,1);
    plot(ts, data(:,1), 'b');
    hold on;
    subplot(3,1,2);
    plot(ts, rad2deg(data(:,2)), 'b');
    hold on;
    subplot(3,1,3);
    plot(ts, rad2deg(data(:,3)), 'b');
    hold on;

    % Plot ellipse error
    lon = rad2deg(sim.sc_pos_vrt.Data(end,2));
    lat = rad2deg(sim.sc_pos_vrt.Data(end,3));
    taem(i,:) = [lon, lat];

    set(0, "CurrentFigure", f_land);
    scatter(lon, lat, 'filled','b');
    hold on;

    % Plot atmospheric density variation
    xdata = sim.rho_rel.Data(:,1) .* 100;
    ydata = (sim.sc_pos_vrt.Data(:,1) - simdb.env.Rp)*1e-3;
    set(0, "CurrentFigure", f_rho);
    plot(xdata, ydata, 'b');
    hold on;

    % Plot AoA & bank angle
    data = sim.sc_ae_ang_gnd.Data;
    set(0, "CurrentFigure", f_alpha);
    plot(ts, rad2deg(data(:,1)), 'b');
    hold on;
    set(0, "CurrentFigure", f_beta);
    plot(ts, rad2deg(data(:,2)), 'b');
    hold on;
    set(0, "CurrentFigure", f_sigma);
    plot(ts, rad2deg(data(:,3)), 'b');
    hold on;

    % Plot nav diff pos
    data = sim.nav_pos_vrt.Data - sim.sc_pos_vrt.Data;
    sz = length(data);
    nav_dR(1:sz,i) = data(:,1);
    nav_dtau(1:sz,i) = rad2deg(data(:,2));
    nav_ddelta(1:sz,i) = rad2deg(data(:,3));
    set(0, "CurrentFigure", f_nav_R);
    subplot(3,1,1);
    plot(ts, data(:,1), 'b');
    hold on;
    subplot(3,1,2);
    plot(ts, rad2deg(data(:,2)), 'b');
    hold on;
    subplot(3,1,3);
    plot(ts, rad2deg(data(:,3)), 'b');
    hold on;

    % Plot nav diff vel
    data = sim.nav_vel_tg.Data - sim.sc_vel_tg.Data;
    sz = length(data);
    nav_dV(1:sz,i) = data(:,1);
    nav_dgamma(1:sz,i) = rad2deg(data(:,2));
    nav_dchi(1:sz,i) = rad2deg(data(:,3));
    set(0, "CurrentFigure", f_nav_V);
    subplot(3,1,1);
    plot(ts, data(:,1), 'b');
    hold on;
    subplot(3,1,2);
    plot(ts, rad2deg(data(:,2)), 'b');
    hold on;
    subplot(3,1,3);
    plot(ts, rad2deg(data(:,3)), 'b');
    hold on;

    % Plot nav diff ae
    data = sim.nav_ae_ang_bod.Data - sim.sc_ae_ang_gnd.Data;
    sz = length(data);
    nav_dalpha(1:sz,i) = rad2deg(data(:,1));
    nav_dbeta(1:sz,i) = rad2deg(data(:,2));
    nav_dsigma(1:sz,i) = rad2deg(data(:,3));
    set(0, "CurrentFigure", f_nav_ae);
    subplot(3,1,1);
    plot(ts, rad2deg(data(:,1)), 'b');
    hold on;
    subplot(3,1,2);
    plot(ts, rad2deg(data(:,2)), 'b');
    hold on;
    subplot(3,1,3);
    plot(ts, rad2deg(data(:,3)), 'b');
    hold on;
end
hold off;

set(0, "CurrentFigure", f_h);
title('Height profile');
xlabel('Time [s]');
ylabel('Height [km]');
grid;

he_mu = mean(height_err,2);
he_std = std(height_err,0,2);
dtau_mu = mean(dtau_err,2);
dtau_std = std(dtau_err,0,2);
ddelta_mu = mean(ddelta_err,2);
ddelta_std = std(ddelta_err,0,2);
set(0, "CurrentFigure", f_dh);
subplot(3,1,1);
hold on;
%plot(max_ts, he_mu, 'yellow', 'LineWidth',3);
plot(max_ts, he_mu - 3*he_std, 'red', 'LineWidth',3);
plot(max_ts, he_mu + 3*he_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta{h} [m]');
grid;
subplot(3,1,2);
hold on;
plot(max_ts, dtau_mu - 3*dtau_std, 'red', 'LineWidth',3);
plot(max_ts, dtau_mu + 3*dtau_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta\tau [deg]');
grid;
subplot(3,1,3);
hold on;
plot(max_ts, ddelta_mu - 3*ddelta_std, 'red', 'LineWidth',3);
plot(max_ts, ddelta_mu + 3*ddelta_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta\delta [deg]');
grid;

set(0, "CurrentFigure", f_V);
title('Velocity profile');
xlabel('Time [s]');
ylabel('Velocity [m/s]');
grid;

Ve_mu = mean(vel_err,2);
Ve_std = std(vel_err,0,2);
dgamma_mu = mean(dgamma_err,2);
dgamma_std = std(dgamma_err,0,2);
dchi_mu = mean(dchi_err,2);
dchi_std = std(dchi_err,0,2);
set(0, "CurrentFigure", f_dV);
subplot(3,1,1);
hold on;
%plot(max_ts, Ve_mu, 'yellow', 'LineWidth',3);
plot(max_ts, Ve_mu - 3*Ve_std, 'red', 'LineWidth',3);
plot(max_ts, Ve_mu + 3*Ve_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta{V} [m/s]');
grid;
subplot(3,1,2);
hold on;
plot(max_ts, dgamma_mu - 3*dgamma_std, 'red', 'LineWidth',3);
plot(max_ts, dgamma_mu + 3*dgamma_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta\gamma [deg]');
grid;
subplot(3,1,3);
hold on;
plot(max_ts, dchi_mu - 3*dchi_std, 'red', 'LineWidth',3);
plot(max_ts, dchi_mu + 3*dchi_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta\chi [deg]');
grid;

set(0, "CurrentFigure", f_land);
hold on;
[a, b] = draw_error_ellipse(remove_outliers(taem,4),3);
hold off;
title('TAEM interface');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
grid;
fprintf('Error ellipse axis: [%f, %f] deg\n', a, b);

set(0, "CurrentFigure", f_rho);
xlabel('Density variation [%]');
ylabel('Height [km]');
ylim([15, 120.2]);
grid;

set(0, "CurrentFigure", f_alpha);
hold on;
plot(sim.AoA_ref.Time, rad2deg(sim.AoA_ref.Data), 'r','LineWidth',3.5);
hold off;
title('Angle of attack profile');
xlabel('Time [s]');
ylabel('\alpha [deg]');
ylim([15,60]);
grid;

set(0, "CurrentFigure", f_beta);
hold on;
plot([0,sim.AoA_ref.Time(end)], [0, 0], 'r','LineWidth',3.5);
hold off;
title('Sideslip angle profile');
xlabel('Time [s]');
ylabel('\beta [deg]');
grid;

set(0, "CurrentFigure", f_sigma);
hold on;
plot(sim.bank_ang_ref.Time, rad2deg(sim.bank_ang_ref.Data), 'r','LineWidth',3.5);
hold off;
title('Bank angle profile');
xlabel('Time [s]');
ylabel('\sigma [deg]');
grid;

% Nav pos delta
dR_mu = mean(nav_dR,2);
dR_std = std(nav_dR,0,2);
dtau_mu = mean(nav_dtau,2);
dtau_std = std(nav_dtau,0,2);
ddelta_mu = mean(nav_ddelta,2);
ddelta_std = std(nav_ddelta,0,2);
set(0, "CurrentFigure", f_nav_R);
subplot(3,1,1);
hold on;
plot(max_ts, dR_mu - 3*dR_std, 'red', 'LineWidth',3);
plot(max_ts, dR_mu + 3*dR_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta{h} [m]');
ylim([-50, 50]);
grid;
subplot(3,1,2);
hold on;
plot(max_ts, dtau_mu - 3*dtau_std, 'red', 'LineWidth',3);
plot(max_ts, dtau_mu + 3*dtau_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta\tau [deg]');
ylim([-5e-4, 5e-4]);
grid;
subplot(3,1,3);
hold on;
plot(max_ts, ddelta_mu - 3*ddelta_std, 'red', 'LineWidth',3);
plot(max_ts, ddelta_mu + 3*ddelta_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta\delta [deg]');
ylim([-5e-4, 5e-4]);
grid;

% Nav vel delta
dV_mu = mean(nav_dV,2);
dV_std = std(nav_dV,0,2);
dgamma_mu = mean(nav_dgamma,2);
dgamma_std = std(nav_dgamma,0,2);
dchi_mu = mean(nav_dchi,2);
dchi_std = std(nav_dchi,0,2);
set(0, "CurrentFigure", f_nav_V);
subplot(3,1,1);
hold on;
plot(max_ts, dV_mu - 3*dV_std, 'red', 'LineWidth',3);
plot(max_ts, dV_mu + 3*dV_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta{V} [m/s]');
ylim([-0.2, 0.2]);
grid;
subplot(3,1,2);
hold on;
plot(max_ts, dgamma_mu - 3*dgamma_std, 'red', 'LineWidth',3);
plot(max_ts, dgamma_mu + 3*dgamma_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta\gamma [deg]');
ylim([-1e-2, 1e-2]);
grid;
subplot(3,1,3);
hold on;
plot(max_ts, dchi_mu - 3*dchi_std, 'red', 'LineWidth',3);
plot(max_ts, dchi_mu + 3*dchi_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta\chi [deg]');
ylim([-5e-3, 5e-3]);
grid;

% Nav ae delta
dalpha_mu = mean(nav_dalpha,2);
dalpha_std = std(nav_dalpha,0,2);
dbeta_mu = mean(nav_dbeta,2);
dbeta_std = std(nav_dbeta,0,2);
dsigma_mu = mean(nav_dsigma,2);
dsigma_std = std(nav_dsigma,0,2);
set(0, "CurrentFigure", f_nav_ae);
subplot(3,1,1);
hold on;
plot(max_ts, dalpha_mu - 3*dalpha_std, 'red', 'LineWidth',3);
plot(max_ts, dalpha_mu + 3*dalpha_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta\alpha [deg]');
%ylim([-0.2, 0.2]);
grid;
subplot(3,1,2);
hold on;
plot(max_ts, dbeta_mu - 3*dbeta_std, 'red', 'LineWidth',3);
plot(max_ts, dbeta_mu + 3*dbeta_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta\beta [deg]');
%ylim([-1e-2, 1e-2]);
grid;
subplot(3,1,3);
hold on;
plot(max_ts, dsigma_mu - 3*dsigma_std, 'red', 'LineWidth',3);
plot(max_ts, dsigma_mu + 3*dsigma_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('\Delta\sigma [deg]');
%ylim([-5e-3, 5e-3]);
grid;