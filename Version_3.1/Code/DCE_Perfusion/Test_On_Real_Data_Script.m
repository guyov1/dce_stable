% Close all open figures
close all;
% Clear everything except breakpoints
current_break_points = dbstatus;
save('myBreakpoints.mat', 'current_break_points');
clear all;
load('myBreakpoints.mat');
delete('myBreakpoints.mat');
dbstop(current_break_points);

% Add current directory and all it's subdirectories to Matlab's path
addpath(genpath('./'));

Verbosity = 'None';

% Initiate simulation struct
Sim_Struct = struct;

% Set simulation parameters
Sim_Struct = Simulation_Set_Params(Sim_Struct, Verbosity);

% Set parallel processing if needed
Set_Parallel_Processing(Sim_Struct, Verbosity);

% Read input MRI data
[Subject_name, Subject_Path, WM_mask_absolute_path, Art_Mask, Vein_Mask, After_CTC_mat, DCECoregP, Brain_Extract_path] = ReadRealData();

%% ---------------------------------------------------------------------------------------------------------------------------------------
PefusionOutput = [Subject_Path filesep 'Run_Output'];
Run_DCE_Perfusion( Subject_name, Subject_Path, Sim_Struct, PefusionOutput, Verbosity );