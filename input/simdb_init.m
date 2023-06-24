%% Construct simulink database
simdb = struct();
simdb = loadInputData(simdb);

%% Process every file inside inputs
function out = loadInputData(simdb)
    noiseFlag = 1;

    out = simdb;
    currdir = pwd;
    files = dir(fullfile(currdir, '**\*.m'));
    
    for i=1:length(files)
        if strcmp(files(i).name, 'simdb_init.m')
            continue
        end
        
        rel_path = erase(files(i).folder, currdir);
        rel_path = rel_path(2:end);
        substr = split(rel_path, '\');
    
        fcn_name = erase(files(i).name, '.m');
        if isempty(substr{1})
            fields = {fcn_name};
        else
            fields = [substr, fcn_name];
        end
        fields = unique(fields, 'stable');
    
        cd(files(i).folder);
        out = setfield(out, fields{:}, eval(fcn_name + "(" + noiseFlag + ")"));
    end
    
    cd(currdir);
end