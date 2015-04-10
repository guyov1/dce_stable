% Pick simulation verbosity (displayed messages indicating simulation progress)
% Optiona : 'None', 'Normal', 'Full'
Verbosity = 'Normal';

% Cleaning before simulation starts
Simulation_Clear;

% Initiate simulation struct
Sim_Struct = struct;

% Set simulation parameters
Sim_Struct = Simulation_Set_Params(Sim_Struct, Verbosity);

% Initiate simulation vars/structs etc.
Sim_Struct = Simulation_Init(Sim_Struct, Verbosity);

% Set parallel processing if needed
Set_Parallel_Processing(Sim_Struct, Verbosity);

if Sim_Struct.RealData_Flag
   error('-E- Running simulation but RealData_Flag is set to TRUE...'); 
end

%% AIFs creation for all iterations
Sim_Struct = Create_AIFs(Sim_Struct, Verbosity);

%% Kernel creation
Sim_Struct = Create_Kernels(Sim_Struct, Verbosity);

%% Driving the differential equation (for a single parameter value)
if (Sim_Struct.Drive_Diff_Eq)
    Sim_Struct = Driving_Differential_Eq(Sim_Struct, Verbosity);
end

%% Filter AIF through kernels
Sim_Struct = Filter_AIF(Sim_Struct, Verbosity);

%% Using Murase to estimate kep, Vp, Ktrans according to Tofts model
if (Sim_Struct.Iterate_Murase_Tofts)
    Sim_Struct = Murase_Tofts_Comparison(Sim_Struct, Verbosity);
end

%% Simulation main body

% Initiate iteration loop vars
num_iterations                  = Sim_Struct.num_iterations;
num_averages                    = Sim_Struct.num_averages;
num_time_stamps                 = Sim_Struct.num_time_stamps;
idx_fig                         = Sim_Struct.idx_fig;
Correct_estimation_due_to_delay = Sim_Struct.Correct_estimation_due_to_delay;
Ignore_Gaussian_Calculation     = Sim_Struct.Ignore_Gaussian_Calculation;
Check_Sourbron_Estimate         = Sim_Struct.Check_Sourbron_Estimate;
num_results_parameters          = Sim_Struct.num_results_parameters;
RealData_Flag                   = Sim_Struct.RealData_Flag;

% Replicate simulation struct and figure index for parallel processing
Sim_Struct_Replicated           = repmat(Sim_Struct,1,num_iterations);
idx_fig_Rep                     = repmat(idx_fig,1,num_iterations);

% Initiate results to zeros
results         = zeros(num_results_parameters,num_iterations);

tic;
if ~strcmp(Verbosity,'None')
    display('-I- STARTED MAIN LOOP!');
end

if (num_iterations == 1 || Sim_Struct.FORCE_SERIAL || Sim_Struct.FORCE_MAIN_LOOP_SERIAL)
    [ results, Sim_Struct_Replicated, idx_fig_Rep] = Simulation_Serial( Sim_Struct_Replicated, idx_fig_Rep, results, num_iterations, num_averages, Verbosity, RealData_Flag); 
else
    [ results, Sim_Struct_Replicated, idx_fig_Rep] = Simulation_Parallel( Sim_Struct_Replicated, idx_fig_Rep, results, num_iterations, num_averages, Verbosity, RealData_Flag); 
end

time_finish = toc;
display(sprintf('-I- Main Loop took %.2f seconds to finish...',time_finish));

% Put the local loop variables in struct
Sim_Struct.results = results;
% Update back the idx figure var
idx_fig            = idx_fig_Rep(1);

%% Results plotting and saving

% Plot simulation results
[ idx_fig ] = Plot_Simulation_Results( Sim_Struct,Verbosity, idx_fig);

% Create PDF Report
Local_Path = [pwd filesep 'Run_Output'];
Log_Path   = [pwd filesep 'Run_Output' filesep 'Log.mat'];
%MakeReport_func(Local_Path, LogFN);
MakeReport_func(Local_Path, Log_Path);

% In case we used only 1 iteration, close all the open figures (we have PDF instead)
if (Sim_Struct.num_iterations == 1)
    close all;
end
close all;
