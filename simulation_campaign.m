N = 30;

if not(isfolder(fullfile(pwd, 'output')))
    mkdir('output')
end

currdate = datetime('now','Format','ddMMyy_HHmmss');
output_folder = 'sim-' + string(currdate);
mkdir(fullfile('output', output_folder));

for i=1:N
    % Load and run simulation i/N
    fprintf('Running simulation case %i/%i...\n', i, N);
    
    run('load_simulation');
    out = sim('simulator',...
              'SimulationMode','accelerator');
    
    % Save input and results to corresponding output folder
    sim_dir = fullfile('output', output_folder, "sim_"+i);
    mkdir(sim_dir);
    save(fullfile(sim_dir, 'input.mat'), 'simdb');
    save(fullfile(sim_dir, 'results.mat'), 'out');

    clear('simdb', 'out');
end