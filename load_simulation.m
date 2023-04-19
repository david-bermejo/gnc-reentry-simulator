%% Initialization, add relevant folders to MATLAB path
addpath("functions\");
addpath("model\lib\");
addpath("model\spacecraft\aerodynamics\");

%% Process every file inside inputs
run input/simdb_init.m