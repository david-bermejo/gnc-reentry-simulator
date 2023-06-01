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
f_h     = figure(1);
f_dh    = figure(2);
f_V     = figure(3);
f_dV    = figure(4);
f_land  = figure(5);
f_rho   = figure(6);

%% Preprocessing
max_sz = 1;
for i=1:N
    sim = results.(['sim',num2str(i)]).out;
    simdb = results.(['input',num2str(i)]).simdb;
    ts = sim.sc_pos_vrt.Time;
    max_sz = max(max_sz, length(ts));
end

%% Create data arrays
taem = zeros(N,2);
height_err = zeros(max_sz,N);
vel_err = zeros(max_sz,N);
max_ts = 0:0.05:(max_sz-1)/20;

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
    data = sim.sc_pos_vrt.Data(:,1) - sim.y_ref.Data(:,1);
    sz = length(data);
    height_err(1:sz,i) = data;
    set(0, "CurrentFigure", f_dh);
    plot(ts, data, 'b');
    hold on;

    % Plot Velocity
    ts = sim.sc_vel_tg.Time;
    data = sim.sc_vel_tg.Data;
    set(0, "CurrentFigure", f_V);
    plot(ts, data(:,1), 'b');
    hold on;

    % Plot Velocity error
    data = sim.sc_vel_tg.Data(:,1) - sim.y_ref.Data(:,4);
    sz = length(data);
    vel_err(1:sz,i) = data;
    set(0, "CurrentFigure", f_dV);
    plot(ts, data, 'b');
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
end
hold off;

set(0, "CurrentFigure", f_h);
xlabel('Time [s]');
ylabel('Height [km]');
grid;

he_mu = mean(height_err,2);
he_std = std(height_err,0,2);
set(0, "CurrentFigure", f_dh);
hold on;
plot(max_ts, he_mu, 'black', 'LineWidth',3);
plot(max_ts, he_mu - 3*he_std, 'red', 'LineWidth',3);
plot(max_ts, he_mu + 3*he_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('Height error [m]');
grid;

set(0, "CurrentFigure", f_V);
xlabel('Time [s]');
ylabel('Velocity [m/s]');
grid;

Ve_mu = mean(vel_err,2);
Ve_std = std(vel_err,0,2);
set(0, "CurrentFigure", f_dV);
hold on;
plot(max_ts, Ve_mu, 'black', 'LineWidth',3);
plot(max_ts, Ve_mu - 3*Ve_std, 'red', 'LineWidth',3);
plot(max_ts, Ve_mu + 3*Ve_std, 'red', 'LineWidth',3);
hold off;
xlabel('Time [s]');
ylabel('Velocity error [m/s]');
grid;

set(0, "CurrentFigure", f_land);
hold on;
[a, b] = draw_error_ellipse(remove_outliers(taem,4),3);
hold off;
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
grid;
fprintf('Error ellipse axis: [%f, %f] deg\n', a, b);

set(0, "CurrentFigure", f_rho);
xlabel('Density variation [%]');
ylabel('Height [km]');
grid;