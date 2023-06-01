N_workers = 6;
N         = 20*N_workers;
M         = N/N_workers;
T         = zeros(1,M);

if not(isfolder(fullfile(pwd, 'output')))
    mkdir('output')
end

currdate = datetime('now','Format','yyMMdd_HHmmss');
output_folder = 'sim-' + string(currdate);
mkdir(fullfile('output', output_folder));
open_system('model/simulator.slx');
for i=1:M
    tic;

    % Load and run parallel simulation i/M
    fprintf('Running simulation case %i/%i...\n', (i-1)*N_workers+1, N);
    
    simdb_arr = cell(1,N_workers);
    for j=1:N_workers
        run('load_simulation');
        in(j) = Simulink.SimulationInput('simulator');
        in(j) = in(j).setModelParameter('SimulationMode','accelerator');
        in(j) = in(j).setVariable('simdb',simdb);
        simdb_arr{j} = simdb;
        clear('simdb');
    end
    out_arr = parsim(in, "UseFastRestart","on");
    
    % Save input and results to corresponding output folder
    for j=1:N_workers
        idx = (i-1)*N_workers+j; 
        sim_dir = fullfile('output', output_folder, "sim_"+idx);
        mkdir(sim_dir);

        simdb = simdb_arr{j};
        out = out_arr(j);
        save(fullfile(sim_dir, 'input.mat'), 'simdb');
        save(fullfile(sim_dir, 'results.mat'), 'out');
        clear('simdb', 'out');
    end
    clear('simdb_arr', 'out_arr');

    T(i) = toc;

    T_hat = mean(T(1:i));
    fprintf('Estimated remaining time: %f min\n\n', T_hat*(M-i)/60);
end