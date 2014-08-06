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

% Add needed functions to parallel workers
if (Sim_Struct.num_iterations > 1)
    
    if strcmp(Verbosity,'Full')
        display('-I- Starting Parallel Processing Setting...');
    end
    
    myPool  = gcp;
    myFiles = {'Summarize_Iteration.m', 'Estimate_ht_Wiener.m', 'Simulation.m', 'Print2Pdf.m', 'gprint.m','AddToLog.m'};
    addAttachedFiles(myPool, myFiles);
    
    if strcmp(Verbosity,'Full')
        display('-I- Finished Parallel Processing Setting...');
    end
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

% Replicate simulation struct and figure index for parallel processing
Sim_Struct_Replicated           = repmat(Sim_Struct,1,num_iterations);
idx_fig_Rep                     = repmat(idx_fig,1,num_iterations);

% Initiate results to zeros
results         = zeros(num_results_parameters,num_iterations);

if strcmp(Verbosity,'Full')
    display('-I- STARTED MAIN LOOP!');
end

if (num_iterations == 1)
    Simulation_Serial;
else
    parfor iter_num = 1 : num_iterations
    %for iter_num = 1 : num_iterations
        
        % Calculate iteration time
        tic;
        
        % Initiate local loop vectors
        num_averages                             = Sim_Struct_Replicated(iter_num).num_averages;
        est_sigma_noise_vec                      = zeros(num_averages,1);
        est_t_d_noise_vec                        = zeros(num_averages,1);
        est_amp_noise_vec                        = zeros(num_averages,1);
        est_F_noise_vec                          = zeros(num_averages,1);
        est_Delay_sec_noise_vec                  = zeros(num_averages,1);
        est_Delay_sec_using_Gaussian_noise_vec   = zeros(num_averages,1);
        est_Ki_Patlak_noise_vec                  = zeros(num_averages,1);
        est_Ki_Two_Comp_noise_vec                = zeros(num_averages,1);
        est_Vb_Patlak_noise_vec                  = zeros(num_averages,1);
        est_E_noise_vec                          = zeros(num_averages,1);
        est_PS_noise_vec                         = zeros(num_averages,1);
        est_Vb_Two_Comp_noise_vec                = zeros(num_averages,1);
        est_Vd_noise_vec                         = zeros(num_averages,1);
        est_Vd_normal_tis_noise_vec              = zeros(num_averages,1);
        est_MTT_noise_vec                        = zeros(num_averages,1);
        est_MTT_normal_tis_noise_vec             = zeros(num_averages,1);
        
        est_F_Sourbron_noise_vec                 = zeros(num_averages,1);
        est_Ki_Sourbron_Two_Comp_noise_vec       = zeros(num_averages,1);
        est_Vb_Sourbron_Two_Comp_noise_vec       = zeros(num_averages,1);
        est_Ve_Sourbron_Two_Comp_noise_vec       = zeros(num_averages,1);
        
        %parfor avg_num = 1 : num_averages
        for avg_num = 1 : num_averages
            
            %% Estimating h(t) by Wiener filter
            [est_gauss_filter_Wiener_noise, est_larss_filter_Wiener_noise,idx_fig_Rep(iter_num) ] = Estimate_ht_Wiener(Sim_Struct_Replicated(iter_num), Verbosity, iter_num, avg_num, idx_fig_Rep(iter_num));
            %%  Estimating h(t) by regularization methods
            [ht_Struct, idx_fig_Rep(iter_num)] = Estimate_ht_Regularization(Sim_Struct_Replicated(iter_num), Verbosity, iter_num, avg_num, idx_fig_Rep(iter_num));
            
            % Correct h(t) estimation if it seems we have delay in AIF
            if Correct_estimation_due_to_delay
                [Sim_Struct_Replicated(iter_num).est_delay_by_AIF_correct, Sim_Struct_Replicated(iter_num).Sim_AIF_with_noise_Regul_shifted,ht_Struct.b_PCA_larss_result_2nd_deriv, idx_fig_Rep(iter_num)] = ...
                    AIF_Delay_Correct(Sim_Struct_Replicated(iter_num), ht_Struct, Verbosity, iter_num, avg_num, idx_fig_Rep(iter_num));
            end
            
            %% Gaussian parameters estimation
            if (~Ignore_Gaussian_Calculation)
                [est_sigma_noise , est_t_d_noise, est_amp_noise,  idx_fig_Rep(iter_num) ] = Estimate_Gauss_Params(Sim_Struct_Replicated(iter_num),ht_Struct, est_gauss_filter_Wiener_noise, Verbosity, iter_num, avg_num, idx_fig_Rep(iter_num));
            else
                est_sigma_noise = NaN;
                est_t_d_noise   = NaN;
                est_amp_noise   = NaN;
            end
            
            %% Larsson parameters estimation
            [Larss_Struct, idx_fig_Rep(iter_num)] = Estimate_Larss_Params(Sim_Struct_Replicated(iter_num), ht_Struct, est_t_d_noise, est_larss_filter_Wiener_noise, Verbosity, iter_num, avg_num, idx_fig_Rep(iter_num));
            
            % Put Larrson Parameters in vectors (according to average_iteration)
            Sim_Struct_Replicated(iter_num).est_F_noise_vec(avg_num)                        = Larss_Struct.est_F_noise;
            Sim_Struct_Replicated(iter_num).est_Delay_sec_noise_vec(avg_num)                = Larss_Struct.est_Delay_sec_noise;
            Sim_Struct_Replicated(iter_num).est_Delay_sec_using_Gaussian_noise_vec(avg_num) = Larss_Struct.est_Delay_sec_using_Gaussian_noise;
            Sim_Struct_Replicated(iter_num).est_Ki_Patlak_noise_vec(avg_num)                = Larss_Struct.est_Ki_Patlak_noise;
            Sim_Struct_Replicated(iter_num).est_Ki_Two_Comp_noise_vec(avg_num)              = Larss_Struct.est_Ki_Two_Comp_noise;
            Sim_Struct_Replicated(iter_num).est_Vb_Patlak_noise_vec(avg_num)                = Larss_Struct.est_Vb_Patlak_noise;
            Sim_Struct_Replicated(iter_num).est_E_noise_vec(avg_num)                        = Larss_Struct.est_E_noise;
            Sim_Struct_Replicated(iter_num).est_PS_noise_vec(avg_num)                       = Larss_Struct.est_PS_noise;
            Sim_Struct_Replicated(iter_num).est_Vb_Two_Comp_noise_vec(avg_num)              = Larss_Struct.est_Vb_Two_Comp_noise;
            Sim_Struct_Replicated(iter_num).est_Vd_noise_vec(avg_num)                       = Larss_Struct.est_Vd_noise;
            Sim_Struct_Replicated(iter_num).est_Vd_normal_tis_noise_vec(avg_num)            = Larss_Struct.est_Vd_normal_tis_noise;
            Sim_Struct_Replicated(iter_num).est_MTT_noise_vec(avg_num)                      = Larss_Struct.est_MTT_noise;
            Sim_Struct_Replicated(iter_num).est_MTT_normal_tis_noise_vec(avg_num)           = Larss_Struct.est_MTT_normal_tis_noise;
            
            Sim_Struct_Replicated(iter_num).est_sigma_noise_vec(avg_num)                    = est_sigma_noise;
            Sim_Struct_Replicated(iter_num).est_t_d_noise_vec(avg_num)                      = est_t_d_noise;
            Sim_Struct_Replicated(iter_num).est_amp_noise_vec(avg_num)                      = est_amp_noise;
            
            if (Check_Sourbron_Estimate)
                Sim_Struct_Replicated(iter_num).est_F_Sourbron_noise_vec(avg_num)               = Larss_Struct.est_F_Two_Comp_Sourbron_noise;
                Sim_Struct_Replicated(iter_num).est_Ki_Sourbron_Two_Comp_noise_vec(avg_num)     = Larss_Struct.est_Ki_Two_Comp_Sourbron_noise;
                Sim_Struct_Replicated(iter_num).est_Vb_Sourbron_Two_Comp_noise_vec(avg_num)     = Larss_Struct.est_Vb_Two_Comp_Sourbron_noise;
                Sim_Struct_Replicated(iter_num).est_Ve_Sourbron_Two_Comp_noise_vec(avg_num)     = Larss_Struct.Ve_Two_Comp_Sourbron_est;
            end
            
        end % end of number of averages
        
        % Summarize iteration results and put in relevant output matrices
        results(:,iter_num) = Summarize_Iteration(Sim_Struct_Replicated(iter_num), Verbosity, iter_num, avg_num);
        
    end
end


% Put the local loop variables in struct
Sim_Struct.results       = results;
% Update back the idx figure var
idx_fig = idx_fig_Rep(1);

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
