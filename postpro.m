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
    set(0, "CurrentFigure", f_dV);
    plot(ts, data, 'b');
    hold on;
end
hold off;

set(0, "CurrentFigure", f_h);
xlabel('Time [s]');
ylabel('Height [km]');
grid;

set(0, "CurrentFigure", f_dh);
xlabel('Time [s]');
ylabel('Height error [m]');
grid;

set(0, "CurrentFigure", f_V);
xlabel('Time [s]');
ylabel('Velocity [m/s]');
grid;

set(0, "CurrentFigure", f_dV);
xlabel('Time [s]');
ylabel('Velocity error [m/s]');
grid;