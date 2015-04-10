function [ ...
    Flow_vec_with_Delay, Flow_vec_no_Delay, Delay_sec_by_Max_Val_with_Delay, Delay_sec_by_Max_Val_no_Delay, ...
    t_delay_single_gauss_sec_vec, sigma_seconds_single_gauss_vec, Amp_single_gauss_vec, ...
    fitted_larsson_with_Delay, fitted_larsson_no_Delay, fitted_larsson_with_Delay_High_F, fitted_larsson_no_Delay_High_F, ...
    fitted_larsson_with_Delay_no_E, fitted_larsson_no_Delay_no_E, fitted_larsson_with_Delay_no_E_High_F, fitted_larsson_no_Delay_no_E_High_F, ...
    fitted_gaussian, fitted_double_gaussian, double_gaussian_param_vec, ...
    Ktrans_with_Delay_vec, Ktrans_no_Delay_vec, ...
    E_with_Delay_vec, E_no_Delay_vec, Ktrans_with_Delay_High_F_vec, Ktrans_no_Delay_High_F_vec, ...
    Vb_with_Delay_vec, Vb_no_Delay_vec, Vb_with_Delay_High_F_vec, Vb_no_Delay_High_F_vec, ...
    Vb_with_Delay_no_E_vec, Vb_no_Delay_no_E_vec, Vb_with_Delay_no_E_High_F_vec, Vb_no_Delay_no_E_High_F_vec, ...
    Ve_with_Delay_vec, Ve_no_Delay_vec, Ve_with_Delay_High_F_vec, Ve_no_Delay_High_F_vec, MTT_with_Delay_vec, MTT_no_Delay_vec, ...
    Ktrans_with_Delay_Patlak_vec, Ktrans_no_Delay_Patlak_vec, Vb_with_Delay_Patlak_vec, Vb_no_Delay_Patlak_vec, ...
    MTT_with_Delay_Patlak_vec, MTT_no_Delay_Patlak_vec ] ...
    = Parallel_Params_Est_Real_Data(Sim_Struct,Est_ht_with_Delay,Est_ht_no_Delay,Ct,AIF_delay_corrected,AIF_no_Delay,idx_fig)

display('--------------------------------------------------------');
display('-I- Starting non linear parameters estimation...');
display('--------------------------------------------------------');

time_vec_minutes_T                    = double(Sim_Struct.time_vec_minutes');
Est_ht_with_Delay_T                   = Est_ht_with_Delay';
Est_ht_no_Delay_T                     = Est_ht_no_Delay';
algorithm_options                     = Sim_Struct.algorithm_options;
LowerBound_Gauss                      = Sim_Struct.LowerBound_Gauss;
UpperBound_Gauss                      = Sim_Struct.UpperBound_Gauss;
LowerBound_Larsson                    = Sim_Struct.LowerBound_Larsson;
UpperBound_Larsson                    = Sim_Struct.UpperBound_Larsson;
Correct_estimation_due_to_delay       = Sim_Struct.Correct_estimation_due_to_delay;
Use_Model_Selection                   = Sim_Struct.Use_Model_Selection;
Ignore_Delay_Model_Selection          = Sim_Struct.Ignore_Delay_Model_Selection;

time_vec_minutes                      = Sim_Struct.time_vec_minutes;
Patlak_Est_Type                       = Sim_Struct.Patlak_Est_Type;
Vb_low                                = Sim_Struct.Vb_low;
RealData_Flag                         = Sim_Struct.RealData_Flag;
USE_ONE_GAUSSIAN                      = Sim_Struct.USE_ONE_GAUSSIAN;
USE_DOUBLE_GAUSSIAN                   = Sim_Struct.USE_DOUBLE_GAUSSIAN;
num_voxels                            = Sim_Struct.num_voxels;
num_time_stamps                       = Sim_Struct.num_time_stamps;
Adjusted_Larsson_Model                = Sim_Struct.Adjusted_Larsson_Model;
min_interval                          = Sim_Struct.min_interval;
init_Ve_guess                         = Sim_Struct.init_Ve_guess;

% Gaussian parameters
t_delay_single_gauss_min_vec          = zeros(1,num_voxels);
t_delay_single_gauss_sec_vec          = zeros(1,num_voxels);
est_var_vec                           = zeros(1,num_voxels);
sigma_in_min_vec                      = zeros(1,num_voxels);
sigma_seconds_single_gauss_vec        = zeros(1,num_voxels);
Amp_single_gauss_vec                  = zeros(1,num_voxels);
double_gaussian_param_vec             = zeros(6,num_voxels);

fitted_larsson_with_Delay             = zeros(num_voxels,num_time_stamps);
fitted_larsson_no_Delay               = zeros(num_voxels,num_time_stamps);
fitted_larsson_with_Delay_High_F      = zeros(num_voxels,num_time_stamps);
fitted_larsson_no_Delay_High_F        = zeros(num_voxels,num_time_stamps);
fitted_larsson_with_Delay_no_E        = zeros(num_voxels,num_time_stamps);
fitted_larsson_no_Delay_no_E          = zeros(num_voxels,num_time_stamps);
fitted_larsson_with_Delay_no_E_High_F = zeros(num_voxels,num_time_stamps);
fitted_larsson_no_Delay_no_E_High_F   = zeros(num_voxels,num_time_stamps);

fitted_gaussian                       = zeros(num_voxels,num_time_stamps);
fitted_double_gaussian                = zeros(num_voxels,num_time_stamps);
% Larsson parameters
Flow_vec_with_Delay                   = zeros(1,num_voxels);
Flow_vec_no_Delay                     = zeros(1,num_voxels);

Delay_sec_by_Max_Val_with_Delay       = zeros(1,num_voxels);
Delay_sec_by_Max_Val_no_Delay         = zeros(1,num_voxels);

Ktrans_with_Delay_vec                 = zeros(1,num_voxels);
Ktrans_no_Delay_vec                   = zeros(1,num_voxels);
Ktrans_with_Delay_Patlak_vec          = zeros(1,num_voxels);
Ktrans_no_Delay_Patlak_vec            = zeros(1,num_voxels);

E_with_Delay_vec                      = zeros(1,num_voxels);
E_no_Delay_vec                        = zeros(1,num_voxels);
Ktrans_with_Delay_High_F_vec          = zeros(1,num_voxels);
Ktrans_no_Delay_High_F_vec            = zeros(1,num_voxels);

Vb_with_Delay_vec                     = zeros(1,num_voxels);
Vb_no_Delay_vec                       = zeros(1,num_voxels);
Vb_with_Delay_High_F_vec              = zeros(1,num_voxels);
Vb_no_Delay_High_F_vec                = zeros(1,num_voxels);
Vb_with_Delay_no_E_vec                = zeros(1,num_voxels);
Vb_no_Delay_no_E_vec                  = zeros(1,num_voxels);
Vb_with_Delay_no_E_High_F_vec         = zeros(1,num_voxels);
Vb_no_Delay_no_E_High_F_vec           = zeros(1,num_voxels);
Vb_with_Delay_Patlak_vec              = zeros(1,num_voxels);
Vb_no_Delay_Patlak_vec                = zeros(1,num_voxels);

Ve_with_Delay_vec                     = zeros(1,num_voxels);
Ve_no_Delay_vec                       = zeros(1,num_voxels);
Ve_with_Delay_High_F_vec              = zeros(1,num_voxels);
Ve_no_Delay_High_F_vec                = zeros(1,num_voxels);

MTT_with_Delay_vec                    = zeros(1,num_voxels);
MTT_no_Delay_vec                      = zeros(1,num_voxels);
MTT_with_Delay_Patlak_vec             = zeros(1,num_voxels);
MTT_no_Delay_Patlak_vec               = zeros(1,num_voxels);


parfor j=1:num_voxels
    
    % lsqcurvefit parameters are:
    % analytic function, initial parameters, time vector, data points ,lower
    % and upper bounds and algorithm options
    
    if (USE_ONE_GAUSSIAN)
        %[est_params, residue_norm, residual, exitflag,algo_info] = ...
        [est_params, ~, ~, ~, ~] = ...
            lsqcurvefit(Gaussian_function,init_guess,time_vec_minutes_T,Est_ht_no_Delay_T(:,j),LowerBound_Gauss,UpperBound_Gauss,algorithm_options);
        
        % Put parameters in wanted variables
        t_delay_single_gauss_min_vec(j)   = est_params(1);
        t_delay_single_gauss_sec_vec(j)   = 60 * t_delay_single_gauss_min_vec(j); % Convert delay time to seconds
        est_var_vec(j)                    = est_params(2);
        sigma_in_min_vec(j)               = sqrt(est_var_vec(j));
        sigma_seconds_single_gauss_vec(j) = 60 * sigma_in_min_vec(j); % Convert sigma to seconds
        Amp_single_gauss_vec(j)           = est_params(3);
        
        % Calculate gaussian out of parameters
        fitted_gaussian(j,:)           = Gaussian(time_vec_minutes, t_delay_single_gauss_min_vec(j), sigma_in_min_vec(j)^2, Amp_single_gauss_vec(j));
        
    end
    
    if (USE_DOUBLE_GAUSSIAN)
        %[est_params_2, residue_norm_2, residual_2, exitflag_2,algo_info_2] = ...
        [est_params_2, ~, ~, ~,~] = ...
            lsqcurvefit(Double_Gaussian_function,init_guess_2,time_vec_minutes_T,Est_ht_no_Delay_T(:,j),LowerBound_2,UpperBound_2,algorithm_options);
        
        % Convert sigma to variance and minutes to seconds
        delay_1_seconds = est_params_2(1) * 60;
        sigma_1_seconds = sqrt(est_params_2(2)) * 60;
        amp_1           = est_params_2(3);
        delay_2_seconds = est_params_2(4) * 60;
        sigma_2_seconds = sqrt(est_params_2(5)) * 60;
        amp_2           =  est_params_2(6);
        
        double_gaussian_param_vec(:,j) = [delay_1_seconds sigma_1_seconds amp_1 delay_2_seconds sigma_2_seconds amp_2];
        
        % Calculate double gaussian out of parameters
        fitted_double_gaussian(j,:) = DoubleGaussian(time_vec_minutes, est_params_2(1), est_params_2(2), est_params_2(3) ...
            ,  est_params_2(4), est_params_2(5), est_params_2(6));
        
    end
    
    % Larsson parameters
    Flow_vec_with_Delay(j)            = max(Est_ht_with_Delay_T(:,j));
    Flow_vec_no_Delay(j)              = max(Est_ht_no_Delay_T(:,j));
    
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
    est_F_noise_with_Delay         = Flow_vec_with_Delay(j);
    est_F_noise_no_Delay           = Flow_vec_no_Delay(j);
    Verbosity                      = 'None';
    iter_num                       = 1;
    avg_num                        = 1;
    
    %Use patlak to get initial parameters estimation
    if Correct_estimation_due_to_delay
        [est_Ktrans_Patlak_noise_with_Delay, est_Vb_Patlak_noise_with_Delay ,...
            est_E_Patlak_noise_with_Delay, est_MTT_Patlak_noise_with_Delay, ~] ...
            = Patlak_Estimation(In_Struct,  AIF_delay_corrected(j,:)', Ct(j,:), est_F_noise_with_Delay, Verbosity, iter_num, avg_num, idx_fig);
        
        if ~Ignore_Delay_Model_Selection
            [est_Ktrans_Patlak_noise_no_Delay, est_Vb_Patlak_noise_no_Delay ,...
                est_E_Patlak_noise_no_Delay, est_MTT_Patlak_noise_no_Delay, ~]     ...
                = Patlak_Estimation(In_Struct,  AIF_no_Delay'       , Ct(j,:), est_F_noise_no_Delay  , Verbosity, iter_num, avg_num, idx_fig);
        end
    else
        % With and without delay is the same
        [est_Ktrans_Patlak_noise_no_Delay, est_Vb_Patlak_noise_no_Delay ,...
            est_E_Patlak_noise_no_Delay, est_MTT_Patlak_noise_no_Delay, ~]     ...
            = Patlak_Estimation(In_Struct,  AIF_no_Delay'        , Ct(j,:), est_F_noise_no_Delay  , Verbosity, iter_num, avg_num, idx_fig);
        
        est_Ktrans_Patlak_noise_with_Delay = est_Ktrans_Patlak_noise_no_Delay;
        est_Vb_Patlak_noise_with_Delay     = est_Vb_Patlak_noise_no_Delay;
        est_E_Patlak_noise_with_Delay      = est_E_Patlak_noise_no_Delay;
        est_MTT_Patlak_noise_with_Delay    = est_MTT_Patlak_noise_no_Delay;
    end
    
    %%
    
    % Initial Guess for non-linear curve fitting for Larsson (Vb, E, Ve)
    
    % Full 2CXM Model, with Delay     (5 parameters, Delay, Flow, Ktrans, Vp, Ve)
    Init_Guess_Larsson_with_Delay             = double( [est_Vb_Patlak_noise_with_Delay est_E_Patlak_noise_with_Delay init_Ve_guess] );
    est_F_with_Delay                          = double(Flow_vec_with_Delay(j));
    % Full 2CXM Model, no Delay       (4 parameters, Flow, Ktrans, Vp, Ve)
    Init_Guess_Larsson_no_Delay               = double( [est_Vb_Patlak_noise_no_Delay est_E_Patlak_noise_no_Delay init_Ve_guess] );
    est_F_no_Delay                            = double(Flow_vec_no_Delay(j));
    % Highly Perfused, F->Inf + Delay (4 parameters, Ktrans, Vp, Ve, Delay) - Same as ETM
    Init_Guess_Larsson_with_Delay_High_F      = double( [est_Vb_Patlak_noise_with_Delay est_Ktrans_Patlak_noise_with_Delay init_Ve_guess] );
    est_F_with_Delay_High_F                   = double(exp(100)); % Infinity
    % Highly Perfused, F->Inf         (3 parameters, Ktrans, Vp, Ve) - Same as ETM
    Init_Guess_Larsson_no_Delay_High_F        = double( [est_Vb_Patlak_noise_no_Delay est_Ktrans_Patlak_noise_no_Delay init_Ve_guess] );
    est_F_no_Delay_High_F                     = double(exp(100)); % Infinity
    % No Indicator Exchange/Highly Vascularized  : Ve->0, E->0 +Delay (3 parameters, Flow, Vp, Delay)
    Init_Guess_Larsson_with_Delay_no_E        = double( est_Vb_Patlak_noise_with_Delay );
    est_F_with_Delay_no_E                     = double(Flow_vec_with_Delay(j));
    % No Indicator Exchange/Highly Vascularized  : Ve->0, E->0 (2 parameters, Flow, Vp)
    Init_Guess_Larsson_no_Delay_no_E          = double( est_Vb_Patlak_noise_no_Delay );
    est_F_no_Delay_no_E                       = double(Flow_vec_no_Delay(j));
    % No Indicator Exchange & Highly Perfused, PS->0 & F->inf + Delay: , (2 parameters, Vp, Delay)
    Init_Guess_Larsson_with_Delay_no_E_High_F = double( est_Vb_Patlak_noise_with_Delay );
    est_F_with_Delay_no_E_High_F              = double(exp(100)); % Infinity
    % No Indicator Exchange & Highly Perfused, PS->0 & F->inf : , (1 parameters, Vp)
    Init_Guess_Larsson_no_Delay_no_E_High_F   = double( est_Vb_Patlak_noise_no_Delay );
    est_F_no_Delay_no_E_High_F                = double(exp(100)); % Infinity
    % 0 parameters will be just the mean value
    
    % The analytic funcational of a Larsson function
    if (Adjusted_Larsson_Model)
        Larsson_function_with_Delay             = @(x,t) Adjusted_Larsson_Filter            ( t, est_F_with_Delay      , x(1), x(2), x(3));
        Larsson_function_no_Delay               = @(x,t) Adjusted_Larsson_Filter            ( t, est_F_no_Delay        , x(1), x(2), x(3));
        Larsson_function_with_Delay_High_F      = @(x,t) Adjusted_Larsson_Filter_High_F     ( t,                         x(1), x(2), x(3));
        Larsson_function_no_Delay_High_F        = @(x,t) Adjusted_Larsson_Filter_High_F     ( t,                         x(1), x(2), x(3));
        Larsson_function_with_Delay_no_E        = @(x,t) Adjusted_Larsson_Filter_no_E       ( t, est_F_with_Delay_no_E , x(1)            );
        Larsson_function_no_Delay_no_E          = @(x,t) Adjusted_Larsson_Filter_no_E       ( t, est_F_no_Delay_no_E   , x(1)            );
        Larsson_function_with_Delay_no_E_High_F = @(x,t) Adjusted_Larsson_Filter_no_E_High_F( t,                         x(1)            );
        Larsson_function_no_Delay_no_E_High_F   = @(x,t) Adjusted_Larsson_Filter_no_E_High_F( t,                         x(1)            );
    else
        Larsson_function_with_Delay             = @(x,t) Larsson_Filter( t, est_F_with_Delay             , x(1), x(2), x(3), Hct);
        Larsson_function_no_Delay               = @(x,t) Larsson_Filter( t, est_F_no_Delay               , x(1), x(2), x(3), Hct);
        Larsson_function_with_Delay_High_F      = @(x,t) Larsson_Filter( t, est_F_with_Delay_High_F      , x(1), x(2), x(3), Hct);
        Larsson_function_no_Delay_High_F        = @(x,t) Larsson_Filter( t, est_F_no_Delay_High_F        , x(1), x(2), x(3), Hct);
        Larsson_function_with_Delay_no_E        = @(x,t) Larsson_Filter( t, est_F_with_Delay             , x(1), 0   , x(3), Hct);
        Larsson_function_no_Delay_no_E          = @(x,t) Larsson_Filter( t, est_F_no_Delay               , x(1), 0   , x(3), Hct);
        Larsson_function_with_Delay_no_E_High_F = @(x,t) Larsson_Filter( t, est_F_with_Delay_no_E_High_F , x(1), 0   , x(3), Hct);
        Larsson_function_no_Delay_no_E_High_F   = @(x,t) Larsson_Filter( t, est_F_no_Delay_no_E_High_F   , x(1), 0   , x(3), Hct);
    end
    
    
    % Calculate parameters only in case F is not 0 (or else there is a
    % problem)
    
    if (est_F_no_Delay~=0)
        
        if Correct_estimation_due_to_delay
            
            % Fitting with delay
            %[est_params_Larsson_with_Delay_noise,residue_norm_Larsson_with_Delay_noise,residual_Larsson_with_Delay_noise,exitflag_Larsson_with_Delay_noise,algo_info_Larsson_with_Delay_noise] = ...
            [est_params_Larsson_with_Delay_noise,~,~,~,~] = ...
                lsqcurvefit(Larsson_function_with_Delay,Init_Guess_Larsson_with_Delay,time_vec_minutes_T,Est_ht_with_Delay_T(:,j)/Flow_vec_with_Delay(j),...
                LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
            
            % Fitting with no delay
            if ~Ignore_Delay_Model_Selection
                
                %[est_params_Larsson_no_Delay_noise,residue_norm_Larsson_no_Delay_noise,residual_Larsson_no_Delay_noise,exitflag_Larsson_no_Delay_noise,algo_info_Larsson_no_Delay_noise] = ...
                [est_params_Larsson_no_Delay_noise,~,~,~,~] = ...
                    lsqcurvefit(Larsson_function_no_Delay,Init_Guess_Larsson_no_Delay,time_vec_minutes_T,Est_ht_no_Delay_T(:,j)/Flow_vec_no_Delay(j),...
                    LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
                
            end
            
        else
            % With and without delay is the same
            
            % Fitting with delay
            %[est_params_Larsson_with_Delay_noise,residue_norm_Larsson_with_Delay_noise,residual_Larsson_with_Delay_noise,exitflag_Larsson_with_Delay_noise,algo_info_Larsson_with_Delay_noise] = ...
            [est_params_Larsson_with_Delay_noise,~,~,~,~] = ...
                lsqcurvefit(Larsson_function_with_Delay,Init_Guess_Larsson_with_Delay,time_vec_minutes_T,Est_ht_with_Delay_T(:,j)/Flow_vec_with_Delay(j),...
                LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
            
            % Fitting with no delay
            est_params_Larsson_no_Delay_noise   = est_params_Larsson_with_Delay_noise;
            %residue_norm_Larsson_no_Delay_noise = residue_norm_Larsson_with_Delay_noise;
            %residual_Larsson_no_Delay_noise     = residual_Larsson_with_Delay_noise;
            %exitflag_Larsson_no_Delay_noise     = exitflag_Larsson_with_Delay_noise;
            %algo_info_Larsson_no_Delay_noise    = algo_info_Larsson_with_Delay_noise;
            
        end
        
        % Assigning two compartment parameters estimation
        Vb_with_Delay_vec(j)         = est_params_Larsson_with_Delay_noise(1);
        Vb_no_Delay_vec(j)           = est_params_Larsson_no_Delay_noise(1);
        E_with_Delay_vec(j)          = est_params_Larsson_with_Delay_noise(2);
        E_no_Delay_vec(j)            = est_params_Larsson_no_Delay_noise(2);
        Ktrans_with_Delay_vec(j)     = est_params_Larsson_with_Delay_noise(2)* Flow_vec_with_Delay(j);
        Ktrans_no_Delay_vec(j)       = est_params_Larsson_no_Delay_noise(2)* Flow_vec_no_Delay(j);
        Ve_with_Delay_vec(j)         = est_params_Larsson_with_Delay_noise(3);
        Ve_no_Delay_vec(j)           = est_params_Larsson_no_Delay_noise(3);
        
        % Estimate MTT
        est_IRF_with_Delay           = Est_ht_with_Delay_T(:,j) / Flow_vec_with_Delay(j);
        est_IRF_no_Delay             = Est_ht_no_Delay_T(:,j) / Flow_vec_no_Delay(j);
        
        est_MTT_noise_with_Delay     = cumtrapz(time_vec_minutes,est_IRF_with_Delay);
        MTT_with_Delay_vec(j)        = est_MTT_noise_with_Delay(end);
        
        est_MTT_noise_no_Delay       = cumtrapz(time_vec_minutes,est_IRF_no_Delay);
        MTT_no_Delay_vec(j)          = est_MTT_noise_no_Delay(end);
        
        if (Adjusted_Larsson_Model)
            fitted_larsson_with_Delay(j,:) = Flow_vec_with_Delay(j)*Adjusted_Larsson_Filter( time_vec_minutes, Flow_vec_with_Delay(j), Vb_with_Delay_vec(j), E_with_Delay_vec(j), Ve_with_Delay_vec(j));
            fitted_larsson_no_Delay(j,:)   = Flow_vec_no_Delay(j)  *Adjusted_Larsson_Filter( time_vec_minutes, Flow_vec_no_Delay(j), Vb_no_Delay_vec(j), E_no_Delay_vec(j), Ve_no_Delay_vec(j));
        else
            fitted_larsson_with_Delay(j,:) = Flow_vec_with_Delay(j)*Larsson_Filter( time_vec_minutes, Flow_vec_with_Delay(j), Vb_with_Delay_vec(j), E_with_Delay_vec(j), Ve_with_Delay_vec(j), Hct);
            fitted_larsson_no_Delay(j,:)   = Flow_vec_no_Delay(j)  *Larsson_Filter( time_vec_minutes, Flow_vec_no_Delay(j), Vb_no_Delay_vec(j), E_no_Delay_vec(j), Ve_no_Delay_vec(j), Hct);
        end
        
        if Use_Model_Selection
            
            % ------------------------------------- Fit with high f + Delay ---------------------------------------------------------------
            %To_Fit_High_F_with_Delay = Est_ht_with_Delay_T(:,j)/Flow_vec_with_Delay(j);
            To_Fit_High_F_with_Delay = Est_ht_with_Delay_T(:,j);
            %[est_params_Larsson_noise_with_Delay_High_F,residue_norm_Larsson_noise_with_Delay_High_F,residual_Larsson_noise_with_Delay_High_F,exitflag_Larsson_noise_with_Delay_High_F,algo_info_Larsson_noise_with_Delay_High_F] = ...
            [est_params_Larsson_noise_with_Delay_High_F,~,~,~,~] = ...
                lsqcurvefit(Larsson_function_with_Delay_High_F,Init_Guess_Larsson_with_Delay_High_F,time_vec_minutes_T,To_Fit_High_F_with_Delay,...
                LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
            
            % Assigning two compartment parameters estimation
            Vb_with_Delay_High_F_vec(j)       = est_params_Larsson_noise_with_Delay_High_F(1);
            Ktrans_with_Delay_High_F_vec(j)        = est_params_Larsson_noise_with_Delay_High_F(2);
            Ve_with_Delay_High_F_vec(j)       = est_params_Larsson_noise_with_Delay_High_F(3);
            
            % ------------------------------------- Fit with high f ---------------------------------------------------------------
            if ~Ignore_Delay_Model_Selection
                
                %To_Fit_High_F_no_Delay = Est_ht_no_Delay_T(:,j)/Flow_vec_no_Delay(j);
                To_Fit_High_F_no_Delay = Est_ht_no_Delay_T(:,j);
                %[est_params_Larsson_noise_no_Delay_High_F,residue_norm_Larsson_noise_no_Delay_High_F,residual_Larsson_noise_no_Delay_High_F,exitflag_Larsson_noise_no_Delay_High_F,algo_info_Larsson_noise_no_Delay_High_F] = ...
                [est_params_Larsson_noise_no_Delay_High_F,~,~,~,~] = ...
                    lsqcurvefit(Larsson_function_no_Delay_High_F,Init_Guess_Larsson_no_Delay_High_F,time_vec_minutes_T,To_Fit_High_F_no_Delay,...
                    LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
                
                % Assigning two compartment parameters estimation
                Vb_no_Delay_High_F_vec(j)       = est_params_Larsson_noise_no_Delay_High_F(1);
                Ktrans_no_Delay_High_F_vec(j)        = est_params_Larsson_noise_no_Delay_High_F(2);
                Ve_no_Delay_High_F_vec(j)       = est_params_Larsson_noise_no_Delay_High_F(3);
                
            end
            % ------------------------------------- Fit with delay and no permeability ---------------------------------------------------------------
            To_Fit_no_E_with_Delay = Est_ht_with_Delay_T(:,j)/Flow_vec_with_Delay(j);
            
            %[est_params_Larsson_noise_with_Delay_no_E,residue_norm_Larsson_noise_with_Delay_no_E,residual_Larsson_noise_with_Delay_no_E,exitflag_Larsson_noise_with_Delay_no_E,algo_info_Larsson_noise_with_Delay_no_E] = ...
            [est_params_Larsson_noise_with_Delay_no_E,~,~,~,~] = ...
                lsqcurvefit(Larsson_function_with_Delay_no_E,Init_Guess_Larsson_with_Delay_no_E,time_vec_minutes_T,To_Fit_no_E_with_Delay,...
                LowerBound_Larsson(1),UpperBound_Larsson(1),algorithm_options);
            
            % Assigning two compartment parameters estimation
            Vb_with_Delay_no_E_vec(j)         = est_params_Larsson_noise_with_Delay_no_E(1);
            
            % ------------------------------------- Fit with no delay and no permeability ---------------------------------------------------------------
            if ~Ignore_Delay_Model_Selection
                
                To_Fit_no_E_no_Delay = Est_ht_no_Delay_T(:,j)/Flow_vec_no_Delay(j);
                
                %[est_params_Larsson_noise_no_Delay_no_E,residue_norm_Larsson_noise_no_Delay_no_E,residual_Larsson_noise_no_Delay_no_E,exitflag_Larsson_noise_no_Delay_no_E,algo_info_Larsson_noise_no_Delay_no_E] = ...
                [est_params_Larsson_noise_no_Delay_no_E,~,~,~,~] = ...
                    lsqcurvefit(Larsson_function_no_Delay_no_E,Init_Guess_Larsson_no_Delay_no_E,time_vec_minutes_T,To_Fit_no_E_no_Delay,...
                    LowerBound_Larsson(1),UpperBound_Larsson(1),algorithm_options);
                
                % Assigning two compartment parameters estimation
                Vb_no_Delay_no_E_vec(j)         = est_params_Larsson_noise_no_Delay_no_E(1);
                
            end
            % ------------------------------------- Fit with high f and no permeability + Delay ---------------------------------------------------------------
            %To_Fit_no_E_High_F_with_Delay = Est_ht_with_Delay_T(:,j)/Flow_vec_with_Delay(j);
            To_Fit_no_E_High_F_with_Delay = Est_ht_with_Delay_T(:,j);
            %[est_params_Larsson_noise_with_Delay_no_E_High_F,residue_norm_Larsson_noise_with_Delay_no_E_High_F,residual_Larsson_noise_with_Delay_no_E_High_F,exitflag_Larsson_noise_with_Delay_no_E_High_F,algo_info_Larsson_noise_with_Delay_no_E_High_F] = ...
            [est_params_Larsson_noise_with_Delay_no_E_High_F,~,~,~,~] = ...
                lsqcurvefit(Larsson_function_with_Delay_no_E_High_F,Init_Guess_Larsson_with_Delay_no_E_High_F,time_vec_minutes_T,To_Fit_no_E_High_F_with_Delay,...
                LowerBound_Larsson(1),UpperBound_Larsson(1),algorithm_options);
            
            % Assigning two compartment parameters estimation
            Vb_with_Delay_no_E_High_F_vec(j)         = est_params_Larsson_noise_with_Delay_no_E_High_F(1);
            
            % ------------------------------------- Fit with high f and no permeability ---------------------------------------------------------------
            %To_Fit_no_E_High_F_no_Delay = Est_ht_no_Delay_T(:,j)/Flow_vec_no_Delay(j);
            To_Fit_no_E_High_F_no_Delay = Est_ht_no_Delay_T(:,j);
            if ~Ignore_Delay_Model_Selection
                %[est_params_Larsson_noise_no_Delay_no_E_High_F,residue_norm_Larsson_noise_no_Delay_no_E_High_F,residual_Larsson_noise_no_Delay_no_E_High_F,exitflag_Larsson_noise_no_Delay_no_E_High_F,algo_info_Larsson_noise_no_Delay_no_E_High_F] = ...
                [est_params_Larsson_noise_no_Delay_no_E_High_F,~,~,~,~] = ...
                    lsqcurvefit(Larsson_function_no_Delay_no_E_High_F,Init_Guess_Larsson_no_Delay_no_E_High_F,time_vec_minutes_T,To_Fit_no_E_High_F_no_Delay,...
                    LowerBound_Larsson(1),UpperBound_Larsson(1),algorithm_options);
            end
            
            % Assigning two compartment parameters estimation
            Vb_no_Delay_no_E_High_F_vec(j)         = est_params_Larsson_noise_no_Delay_no_E_High_F(1);
            
            if (Adjusted_Larsson_Model)
                fitted_larsson_with_Delay_High_F(j, :)      = Adjusted_Larsson_Filter_High_F      ( time_vec_minutes,                          Vb_with_Delay_High_F_vec(j)      , Ktrans_with_Delay_High_F_vec(j) , Ve_with_Delay_High_F_vec(j));
                if ~Ignore_Delay_Model_Selection
                    fitted_larsson_no_Delay_High_F(j, :)        = Adjusted_Larsson_Filter_High_F      ( time_vec_minutes,                          Vb_no_Delay_High_F_vec(j)        , Ktrans_no_Delay_High_F_vec(j)   , Ve_no_Delay_High_F_vec(j));
                end
                fitted_larsson_with_Delay_no_E(j, :)        = Flow_vec_with_Delay(j) * Adjusted_Larsson_Filter_no_E        ( time_vec_minutes, Flow_vec_with_Delay(j) , Vb_with_Delay_no_E_vec(j)        );
                if ~Ignore_Delay_Model_Selection
                    fitted_larsson_no_Delay_no_E(j, :)          = Flow_vec_no_Delay(j)   * Adjusted_Larsson_Filter_no_E        ( time_vec_minutes, Flow_vec_no_Delay(j)   , Vb_no_Delay_no_E_vec(j)          );
                end
                fitted_larsson_with_Delay_no_E_High_F(j, :) = Adjusted_Larsson_Filter_no_E_High_F ( time_vec_minutes,                          Vb_with_Delay_no_E_High_F_vec(j) );
                if ~Ignore_Delay_Model_Selection
                    fitted_larsson_no_Delay_no_E_High_F(j, :)   = Adjusted_Larsson_Filter_no_E_High_F ( time_vec_minutes,                          Vb_no_Delay_no_E_High_F_vec(j)   );
                end
            else
                fitted_larsson_with_Delay_High_F(j, :)      = Larsson_Filter( time_vec_minutes, est_F_with_Delay_High_F      , Vb_with_Delay_High_F_vec(j)      , 0, Ve_with_Delay_High_F_vec(j), Hct);
                fitted_larsson_no_Delay_High_F(j, :)        = Larsson_Filter( time_vec_minutes, est_F_no_Delay_High_F        , Vb_no_Delay_High_F_vec(j)        , 0, Ve_no_Delay_High_F_vec(j)  , Hct);
                fitted_larsson_with_Delay_no_E(j, :)        = Flow_vec_with_Delay(j) * Larsson_Filter( time_vec_minutes, Flow_vec_with_Delay(j)       , Vb_with_Delay_no_E_vec(j)        , 0, 0                          , Hct);
                fitted_larsson_no_Delay_no_E(j, :)          = Flow_vec_no_Delay(j)   * Larsson_Filter( time_vec_minutes, Flow_vec_no_Delay(j)         , Vb_no_Delay_no_E_vec(j)          , 0, 0                          , Hct);
                fitted_larsson_with_Delay_no_E_High_F(j, :) = Larsson_Filter( time_vec_minutes, est_F_with_Delay_no_E_High_F , Vb_with_Delay_no_E_High_F_vec(j) , 0, 0                          , Hct);
                fitted_larsson_no_Delay_no_E_High_F(j, :)   = Larsson_Filter( time_vec_minutes, est_F_no_Delay_no_E_High_F   , Vb_no_Delay_no_E_High_F_vec(j)   , 0, 0                          , Hct);
            end
        end
        
        % Assigning Patlak parameters estimation
        Ktrans_with_Delay_Patlak_vec(j)    = est_Ktrans_Patlak_noise_with_Delay;
        Ktrans_no_Delay_Patlak_vec(j)      = est_Ktrans_Patlak_noise_no_Delay;
        Vb_with_Delay_Patlak_vec(j)        = est_Vb_Patlak_noise_with_Delay;
        Vb_no_Delay_Patlak_vec(j)          = est_Vb_Patlak_noise_no_Delay;
        MTT_with_Delay_Patlak_vec(j)       = est_MTT_Patlak_noise_with_Delay;
        MTT_no_Delay_Patlak_vec(j)         = est_MTT_Patlak_noise_no_Delay;
        
        % Delay of the AIF will be calculated according to the place of the maximum value of F*IRF
        max_index_with_Delay               = find( Est_ht_with_Delay_T(:,j) == Flow_vec_with_Delay(j) );
        max_index_no_Delay                 = find( Est_ht_no_Delay_T(:,j)   == Flow_vec_no_Delay(j) );
        
        % Translate to minutes
        Delay_sec_by_Max_Val_with_Delay(j) = (max_index_with_Delay - 1) * min_interval * 60;
        Delay_sec_by_Max_Val_no_Delay(j)   = (max_index_no_Delay   - 1) * min_interval * 60;
        
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


