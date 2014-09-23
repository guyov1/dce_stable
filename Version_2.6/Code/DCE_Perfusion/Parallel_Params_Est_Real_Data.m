function [ Flow_vec, Delay_sec_by_Max_Val, t_delay_single_gauss_sec_vec, sigma_seconds_single_gauss_vec, ...
           Amp_single_gauss_vec, fitted_gaussian, fitted_double_gaussian, double_gaussian_param_vec, Ktrans_vec, E_vec, Vb_vec, ...
           Ve_vec, MTT_vec, Ktrans_Patlak_vec, Vb_Patlak_vec, MTT_Patlak_vec ] = Parallel_Params_Est_Real_Data(Sim_Struct, Est_ht, Ct, AIF, idx_fig )

display('--------------------------------------------------------');
display('-I- Starting non linear parameters estimation...');
display('--------------------------------------------------------');

time_vec_minutes_T             = double(Sim_Struct.time_vec_minutes');
Est_ht_T                       = Est_ht';
algorithm_options              = Sim_Struct.algorithm_options;
LowerBound_Gauss               = Sim_Struct.LowerBound_Gauss;
UpperBound_Gauss               = Sim_Struct.UpperBound_Gauss;
LowerBound_Larsson             = Sim_Struct.LowerBound_Larsson;
UpperBound_Larsson             = Sim_Struct.UpperBound_Larsson;

time_vec_minutes               = Sim_Struct.time_vec_minutes;
Patlak_Est_Type                = Sim_Struct.Patlak_Est_Type; 
Vb_low                         = Sim_Struct.Vb_low;
RealData_Flag                  = Sim_Struct.RealData_Flag; 
USE_ONE_GAUSSIAN               = Sim_Struct.USE_ONE_GAUSSIAN;
USE_DOUBLE_GAUSSIAN            = Sim_Struct.USE_DOUBLE_GAUSSIAN;
num_voxels                     = Sim_Struct.num_voxels;
num_time_stamps                = Sim_Struct.num_time_stamps;
Adjusted_Larsson_Model         = Sim_Struct.Adjusted_Larsson_Model;
min_interval                   = Sim_Struct.min_interval;
init_Ve_guess                  = Sim_Struct.init_Ve_guess;

% Gaussian parameters
t_delay_single_gauss_min_vec   = zeros(1,num_voxels);
t_delay_single_gauss_sec_vec   = zeros(1,num_voxels);
est_var_vec                    = zeros(1,num_voxels);
sigma_in_min_vec               = zeros(1,num_voxels);
sigma_seconds_single_gauss_vec = zeros(1,num_voxels);
Amp_single_gauss_vec           = zeros(1,num_voxels);
double_gaussian_param_vec      = zeros(6,num_voxels);
fitted_gaussian                = zeros(num_voxels,num_time_stamps);
fitted_double_gaussian         = zeros(num_voxels,num_time_stamps);
% Larsson parameters
Flow_vec                       = zeros(1,num_voxels);
Delay_sec_by_Max_Val           = zeros(1,num_voxels);
Ktrans_vec                     = zeros(1,num_voxels);
E_vec                          = zeros(1,num_voxels);
Vb_vec                         = zeros(1,num_voxels);
Ve_vec                         = zeros(1,num_voxels);
MTT_vec                        = zeros(1,num_voxels);
Ktrans_Patlak_vec              = zeros(1,num_voxels);
Vb_Patlak_vec                  = zeros(1,num_voxels);
MTT_Patlak_vec                 = zeros(1,num_voxels);


parfor j=1:num_voxels
    
    % lsqcurvefit parameters are:
    % analytic function, initial parameters, time vector, data points ,lower
    % and upper bounds and algorithm options
    
    if (USE_ONE_GAUSSIAN)
        [est_params, residue_norm, residual, exitflag,algo_info] = ...
            lsqcurvefit(Gaussian_function,init_guess,time_vec_minutes_T,Est_ht_T(:,j),LowerBound_Gauss,UpperBound_Gauss,algorithm_options);
    else
        est_params                    = zeros(1,3);
    end
      
    if (USE_DOUBLE_GAUSSIAN)
        [est_params_2, residue_norm_2, residual_2, exitflag_2,algo_info_2] = ...
            lsqcurvefit(Double_Gaussian_function,init_guess_2,time_vec_minutes_T,Est_ht_T(:,j),LowerBound_2,UpperBound_2,algorithm_options);
    else
        est_params_2                  = zeros(1,6);
    end
    
    % Put parameters in wanted variables
    t_delay_single_gauss_min_vec(j)   = est_params(1);
    t_delay_single_gauss_sec_vec(j)   = 60 * t_delay_single_gauss_min_vec(j); % Convert delay time to seconds
    est_var_vec(j)                    = est_params(2);
    sigma_in_min_vec(j)               = sqrt(est_var_vec(j));
    sigma_seconds_single_gauss_vec(j) = 60 * sigma_in_min_vec(j); % Convert sigma to seconds
    Amp_single_gauss_vec(j)           = est_params(3);
    
    % Larsson parameters
    Flow_vec(j)                       = max(Est_ht_T(:,j));
    
    %% ----------------------- PATLAK --------------------------------------
    In_Struct                      = struct;
    In_Struct.time_vec_minutes     = time_vec_minutes;
    %In_Struct.Sim_AIF_with_noise   = AIF;
    %In_Struct.Sim_Ct_larss_kernel  = Ct(j,:);
    In_Struct.plot_flag            = false;
    In_Struct.Ktrans               = NaN; % Simulation ground truth values
    In_Struct.Vb_larss             = NaN; % Simulation ground truth values
    In_Struct.Patlak_Est_Type      = Patlak_Est_Type;
    In_Struct.Vb_low               = Vb_low;
    In_Struct.RealData_Flag        = RealData_Flag;
    est_F_noise                    = Flow_vec(j);
    Verbosity                      = 'None';
    iter_num                       = 1;
    avg_num                        = 1;
    
    %Use patlak to get initial parameters estimation
    [est_Ktrans_Patlak_noise, est_Vb_Patlak_noise ,est_E_Patlak_noise, est_MTT_Patlak_noise, ~] = Patlak_Estimation(In_Struct,  AIF(j,:)', Ct(j,:), est_F_noise, Verbosity, iter_num, avg_num, idx_fig);
    
    %%
    
    % Initial Guess for non-linear curve fitting for Larsson (Vb, E, Ve)
    Init_Guess_Larsson = double( [est_Vb_Patlak_noise est_E_Patlak_noise init_Ve_guess] );
    
    est_F = double(Flow_vec(j));
    % The analytic funcational of a Larsson function
    if (Adjusted_Larsson_Model)
        Larsson_function         = @(x,t) Adjusted_Larsson_Filter( t, est_F, x(1), x(2), x(3));
    else
        Larsson_function         = @(x,t) Larsson_Filter( t, est_F, x(1), x(2), x(3), Hct);
    end
    
    % Calculate parameters only in case F is not 0 (or else there is a
    % problem)
    
    if (est_F~=0)
        [est_params_Larsson_noise,residue_norm_Larsson_noise,residual_Larsson_noise,exitflag_Larsson_noise,algo_info_Larsson_noise] = ...
            lsqcurvefit(Larsson_function,Init_Guess_Larsson,time_vec_minutes_T,Est_ht_T(:,j)/Flow_vec(j),...
            LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
        
        % Assigning two compartment parameters estimation
        Vb_vec(j)         = est_params_Larsson_noise(1);
        E_vec(j)          = est_params_Larsson_noise(2);
        Ktrans_vec(j)     = est_params_Larsson_noise(2)* Flow_vec(j);
        Ve_vec(j)         = est_params_Larsson_noise(3);
        % Estimate MTT
        est_IRF           = Est_ht_T(:,j) / Flow_vec(j);
        est_MTT_noise     = cumtrapz(time_vec_minutes,est_IRF);
        MTT_vec(j)        = est_MTT_noise(end);
        
        % Assigning Patlak parameters estimation
        Ktrans_Patlak_vec(j)  = est_Ktrans_Patlak_noise;
        Vb_Patlak_vec(j)      = est_Vb_Patlak_noise;
        MTT_Patlak_vec(j)     = est_MTT_Patlak_noise;
        
        % Delay of the AIF will be calculated according to the place of the maximum value of F*IRF
        max_index        = find( Est_ht_T(:,j) == Flow_vec(j) );
        % Translate to minutes
        Delay_sec_by_Max_Val(j) =  (max_index - 1) * min_interval * 60;
        
        % Convert sigma to variance and minutes to seconds
        delay_1_seconds = est_params_2(1) * 60;
        sigma_1_seconds = sqrt(est_params_2(2)) * 60;
        amp_1           = est_params_2(3);
        delay_2_seconds = est_params_2(4) * 60;
        sigma_2_seconds = sqrt(est_params_2(5)) * 60;
        amp_2           =  est_params_2(6);
        
        double_gaussian_param_vec(:,j) = [delay_1_seconds sigma_1_seconds amp_1 delay_2_seconds sigma_2_seconds amp_2];
        
        % Calculate gaussian out of parameters
        fitted_gaussian(j,:)           = Gaussian(time_vec_minutes, t_delay_single_gauss_min_vec(j), sigma_in_min_vec(j)^2, Amp_single_gauss_vec(j));
        
        % Calculate double gaussian out of parameters
        fitted_double_gaussian(j,:) = DoubleGaussian(time_vec_minutes, est_params_2(1), est_params_2(2), est_params_2(3) ...
            ,  est_params_2(4), est_params_2(5), est_params_2(6));
        
    end
    
    % Report after each 100 voxels
    if ( mod(j,1000) == 0 )
        
        %         duration_of_100 = toc;
        %         tic;
        %
        display(sprintf('Finished lsqcurvefit for %d voxels...',j));
        %
        %         disp(['Calculating 100 iterations took: ' num2str(duration_of_100) ' seconds.']);
        %
        remaining_voxels = num_voxels - j;
        
        fprintf('Number of remaining voxels: %d .\n',remaining_voxels);
        
        %         % Increment every 100 voxels
        %         sub_iter_num = sub_iter_num + 1;
        %
        %         % Update average each iteration (X_j_mean = X_j-1_mean*(N-1/N) + X_j/N)
        %         duration_of_100_average = duration_of_100_average*((sub_iter_num - 1) / sub_iter_num ) + (duration_of_100/sub_iter_num);
        %
        %         WillFinishAround = start_time + ( ((num_voxels - j)/100)*duration_of_100_average ) /24/60/60;
        %         disp(['Will finish around: ' datestr(WillFinishAround)]);
        
        % Finish calculation time
        %par_iter(j) = Par.toc;
    end
    
    
    
    
end

end


