function [] = Set_Parallel_Processing( Sim_Struct, Verbosity )

% Do not display warning regarding attaching files that were attached before 
warning('off','parallel:lang:pool:IgnoringAlreadyAttachedFiles')

tic;

if ~strcmp(Verbosity,'None')
    display('-I- Setting Parallel Processing parameters...');
end

% Enable parallel processing
if ( ~Sim_Struct.FORCE_SERIAL || ~Sim_Struct.FORCE_MAIN_LOOP_SERIAL )
    try
        num_processes = getenv('NUMBER_OF_PROCESSORS');
        % Maximum allowed processes are 12 in matlab 2012
        if (num_processes > 12)
            num_processes = 11;
        end
        
        myCluster = parcluster('local');
        myCluster.NumWorkers = num_processes; % 'Modified' property now TRUE
        saveProfile(myCluster);   % 'local' profile now updated,
        
        if ~strcmp(Verbosity,'None')
            fprintf('\n');
            display('-I- Initiating matlab pool for parallel processing...');
            fprintf('\n');
        end
        
        %matlabpool;
        parpool;
        
        if ~strcmp(Verbosity,'None')
            fprintf('\n');
            display('-I- Finished matlab pool initiation.');
            fprintf('\n');
            
        end
        
    catch error
        
        if ~strcmp(Verbosity,'None')
            fprintf('\n');
            fprintf('\n');
            display('-I- Matlab pool already running!');
            fprintf('\n');
        end
        
    end
end

% Add needed functions to parallel workers
if ( (Sim_Struct.num_iterations > 1) || ~Sim_Struct.FORCE_SERIAL || ~Sim_Struct.FORCE_MAIN_LOOP_SERIAL || Sim_Struct.RealData_Flag )
    
    if strcmp(Verbosity,'Full')
        display('-I- Starting Parallel Processing Setting...');
    end
    
    myPool  = gcp;
    myFiles = {'Summarize_Iteration.m', 'Estimate_ht_Wiener.m', 'Simulation.m', ...
               'Print2Pdf.m', 'gprint.m', 'AddToLog.m', 'Regularized_Sol.m', ...
               'AIF_Parker.m', 'ReScale_AIF.m', 'Choose_Needed_Ht.m', ...
               'Create_B_matrix.m', 'Basis_spline_function.m', 'PCA_basis.m', ...
               'Patlak_Estimation.m', 'AIF_Delay_Correct.m', 'Gaussian.m', 'DoubleGaussian.m'};
    
    addAttachedFiles(myPool, myFiles);
    
    if strcmp(Verbosity,'Full')
        display('-I- Finished Parallel Processing Setting...');
    end
end

time_finish = toc;

if ~strcmp(Verbosity,'None')
    display(sprintf('-I- Setting Parallel Processing paramteres took %.2f seconds to finish...',time_finish));
end


if strcmp(Verbosity,'Full')
    display('-I- Finished setting Parallel Processing parameters...');
end


end

