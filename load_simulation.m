%% Initialization, add relevant folders to MATLAB path
addpath("functions\");
addpath("model\lib\");

%% Process every file inside input folder and load 'simdb' structure
run input/simdb_init.m