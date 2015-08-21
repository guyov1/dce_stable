function [ gaussian_param, fitted_gaussian, double_gaussian_param, fitted_double_gaussian,...
    Flow_with_Delay, Flow_no_Delay, Vb_with_Delay, Vb_no_Delay, E_with_Delay, E_no_Delay,Ktrans_with_Delay, Ktrans_no_Delay, Ve_with_Delay, Ve_no_Delay, MTT_with_Delay, MTT_no_Delay, fitted_larsson_with_Delay, fitted_larsson_no_Delay,...
    Vb_with_Delay_High_F, Ktrans_with_Delay_High_F, Ve_with_Delay_High_F, fitted_larsson_with_Delay_High_F,...
    Vb_no_Delay_High_F, Ktrans_no_Delay_High_F, Ve_no_Delay_High_F, fitted_larsson_no_Delay_High_F,...
    Vb_with_Delay_no_Ve, E_with_Delay_no_Ve, fitted_larsson_with_Delay_no_Ve,...
    Vb_no_Delay_no_Ve, E_no_Delay_no_Ve, fitted_larsson_no_Delay_no_Ve,...
    Vb_with_Delay_no_E, fitted_larsson_with_Delay_no_E,...
    Vb_no_Delay_no_E, fitted_larsson_no_Delay_no_E,...
    Vb_with_Delay_no_E_High_F, fitted_larsson_with_Delay_no_E_High_F,...
    Vb_no_Delay_no_E_High_F, fitted_larsson_no_Delay_no_E_High_F,...
    Ktrans_with_Delay_Patlak, Ktrans_no_Delay_Patlak, Vb_with_Delay_Patlak, Vb_no_Delay_Patlak, MTT_with_Delay_Patlak, MTT_no_Delay_Patlak, Delay_sec_by_Max_Val_with_Delay, Delay_sec_by_Max_Val_no_Delay] = ....
    Params_Est_Real_Data( Sim_Struct, Est_ht_with_Delay_T, Est_ht_no_Delay_T, Ct, AIF_delay_corrected, AIF_no_Delay, idx_fig)
%Params_Est_Real_Data Estimate parameters by curve fitting

% Take needed parameters
time_vec_minutes_T                    = double(Sim_Struct.time_vec_minutes');
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
Adjusted_Larsson_Model                = Sim_Struct.Adjusted_Larsson_Model;
min_interval                          = Sim_Struct.min_interval;
init_Ve_guess                         = Sim_Struct.init_Ve_guess;
Hct                                   = Sim_Struct.Hct_single;


%% Initiate all return parameters with zeros
gaussian_param          = zeros(3, 1);
fitted_gaussian         = zeros(size(time_vec_minutes));
double_gaussian_param   = zeros(6, 1);
fitted_double_gaussian  = zeros(size(time_vec_minutes));

Flow_with_Delay           = 0;
Flow_no_Delay             = 0;
Vb_with_Delay             = 0;
Vb_no_Delay               = 0;
E_with_Delay              = 0;
E_no_Delay_no_Ve          = 0;
E_with_Delay_no_Ve        = 0;
E_no_Delay                = 0;
Ktrans_with_Delay         = 0;
Ktrans_no_Delay           = 0;
Ve_with_Delay             = 0;
Ve_no_Delay               = 0;
MTT_with_Delay            = 0;
MTT_no_Delay              = 0;
fitted_larsson_with_Delay = zeros(size(time_vec_minutes));
fitted_larsson_no_Delay   = zeros(size(time_vec_minutes));

Vb_with_Delay_High_F             = 0;
Vb_with_Delay_no_Ve              = 0;
Ktrans_with_Delay_High_F         = 0;
Ve_with_Delay_High_F             = 0;
fitted_larsson_with_Delay_High_F = zeros(size(time_vec_minutes));
fitted_larsson_with_Delay_no_Ve  = zeros(size(time_vec_minutes));

Vb_no_Delay_High_F             = 0;
Vb_no_Delay_no_Ve              = 0;
Ktrans_no_Delay_High_F         = 0;
Ve_no_Delay_High_F             = 0;
fitted_larsson_no_Delay_High_F = zeros(size(time_vec_minutes));
fitted_larsson_no_Delay_no_Ve  = zeros(size(time_vec_minutes));

Vb_with_Delay_no_E             = 0;
fitted_larsson_with_Delay_no_E = zeros(size(time_vec_minutes));

Vb_no_Delay_no_E               = 0;
fitted_larsson_no_Delay_no_E   = zeros(size(time_vec_minutes));

Vb_with_Delay_no_E_High_F             = 0;
fitted_larsson_with_Delay_no_E_High_F = zeros(size(time_vec_minutes));

Vb_no_Delay_no_E_High_F             = 0;
fitted_larsson_no_Delay_no_E_High_F = zeros(size(time_vec_minutes));

Ktrans_with_Delay_Patlak         = 0;
Ktrans_no_Delay_Patlak           = 0;
Vb_with_Delay_Patlak             = 0;
Vb_no_Delay_Patlak               = 0;
MTT_with_Delay_Patlak            = 0;
MTT_no_Delay_Patlak              = 0;
Delay_sec_by_Max_Val_with_Delay  = 0;
Delay_sec_by_Max_Val_no_Delay    = 0;

% Boundaries for Vp, E
LowerBound_Larsson_no_Ve        = LowerBound_Larsson(1:2);
UpperBound_Larsson_no_Ve        = UpperBound_Larsson(1:2);
% Boundaries for Vp
LowerBound_Larsson_no_E         = LowerBound_Larsson(1);
UpperBound_Larsson_no_E         = UpperBound_Larsson(1);


%% Estimation
% lsqcurvefit parameters are:
% analytic function, initial parameters, time vector, data points ,lower
% and upper bounds and algorithm options

if (USE_ONE_GAUSSIAN)
    %[est_params, residue_norm, residual, exitflag,algo_info] = ...
    [est_params, ~, ~, ~, ~] = lsqcurvefit(Gaussian_function,init_guess,time_vec_minutes_T,Est_ht_with_Delay_T,LowerBound_Gauss,UpperBound_Gauss,algorithm_options);
    
    % Put parameters in wanted variables
    t_delay_single_gauss_min   = est_params(1);
    t_delay_single_gauss_sec   = 60 * t_delay_single_gauss_min; % Convert delay time to seconds
    est_var                    = est_params(2);
    sigma_in_min               = sqrt(est_var);
    sigma_seconds_single_gauss = 60 * sigma_in_min; % Convert sigma to seconds
    Amp_single_gauss           = est_params(3);
    
    gaussian_param             = [t_delay_single_gauss_sec sigma_seconds_single_gauss Amp_single_gauss ];
    
    % Calculate gaussian out of parameters
    fitted_gaussian            = Gaussian(time_vec_minutes, t_delay_single_gauss_min, sigma_in_min^2, Amp_single_gauss);
end

if (USE_DOUBLE_GAUSSIAN)
    %[est_params_2, residue_norm_2, residual_2, exitflag_2,algo_info_2] = ...
    [est_params_2, ~, ~, ~,~] = ...
        lsqcurvefit(Double_Gaussian_function,init_guess_2,time_vec_minutes_T,Est_ht_with_Delay_T,LowerBound_2,UpperBound_2,algorithm_options);
    
    % Convert sigma to variance and minutes to seconds
    delay_1_seconds = est_params_2(1) * 60;
    sigma_1_seconds = sqrt(est_params_2(2)) * 60;
    amp_1           = est_params_2(3);
    delay_2_seconds = est_params_2(4) * 60;
    sigma_2_seconds = sqrt(est_params_2(5)) * 60;
    amp_2           =  est_params_2(6);
    
    double_gaussian_param = [delay_1_seconds sigma_1_seconds amp_1 delay_2_seconds sigma_2_seconds amp_2];
    
    % Calculate double gaussian out of parameters
    fitted_double_gaussian = DoubleGaussian(time_vec_minutes, est_params_2(1), est_params_2(2), est_params_2(3) ...
        ,  est_params_2(4), est_params_2(5), est_params_2(6));
end

% Larsson parameters
Flow_with_Delay            = max(Est_ht_with_Delay_T);
Flow_no_Delay              = max(Est_ht_no_Delay_T);

%% ----------------------- PATLAK --------------------------------------
In_Struct                      = struct;
In_Struct.time_vec_minutes     = time_vec_minutes;
%In_Struct.Sim_AIF_with_noise   = AIF;
%In_Struct.Sim_Ct_larss_kernel  = Ct;
In_Struct.plot_flag            = false;
In_Struct.Ktrans               = NaN; % Simulation ground truth values
In_Struct.Vb_larss             = NaN; % Simulation ground truth values
In_Struct.Patlak_Est_Type      = Patlak_Est_Type;
In_Struct.Vb_low               = Vb_low;
In_Struct.RealData_Flag        = RealData_Flag;
est_F_noise_with_Delay         = Flow_with_Delay;
est_F_noise_no_Delay           = Flow_no_Delay;
Verbosity                      = 'None';
iter_num                       = 1;
avg_num                        = 1;

%Use patlak to get initial parameters estimation
if Correct_estimation_due_to_delay
    [est_Ktrans_Patlak_noise_with_Delay, est_Vb_Patlak_noise_with_Delay ,...
        est_E_Patlak_noise_with_Delay, est_MTT_Patlak_noise_with_Delay, ~] ...
        = Patlak_Estimation(In_Struct,  AIF_delay_corrected, Ct, est_F_noise_with_Delay, Verbosity, iter_num, avg_num, idx_fig);
    
    if ~Ignore_Delay_Model_Selection
        [est_Ktrans_Patlak_noise_no_Delay, est_Vb_Patlak_noise_no_Delay ,...
            est_E_Patlak_noise_no_Delay, est_MTT_Patlak_noise_no_Delay, ~]     ...
            = Patlak_Estimation(In_Struct,  AIF_no_Delay, Ct, est_F_noise_no_Delay  , Verbosity, iter_num, avg_num, idx_fig);
    end
else
    % With and without delay is the same
    [est_Ktrans_Patlak_noise_no_Delay, est_Vb_Patlak_noise_no_Delay ,...
        est_E_Patlak_noise_no_Delay, est_MTT_Patlak_noise_no_Delay, ~]     ...
        = Patlak_Estimation(In_Struct,  AIF_no_Delay, Ct, est_F_noise_no_Delay  , Verbosity, iter_num, avg_num, idx_fig);
    
    est_Ktrans_Patlak_noise_with_Delay = est_Ktrans_Patlak_noise_no_Delay;
    est_Vb_Patlak_noise_with_Delay     = est_Vb_Patlak_noise_no_Delay;
    est_E_Patlak_noise_with_Delay      = est_E_Patlak_noise_no_Delay;
    est_MTT_Patlak_noise_with_Delay    = est_MTT_Patlak_noise_no_Delay;
end

%%

% Initial Guess for non-linear curve fitting for Larsson (Vb, E, Ve)

% Full 2CXM Model, with Delay     (5 parameters, Delay, Flow, Ktrans, Vp, Ve)
Init_Guess_Larsson_with_Delay             = double( [est_Vb_Patlak_noise_with_Delay est_E_Patlak_noise_with_Delay init_Ve_guess] );
est_F_with_Delay                          = double(Flow_with_Delay);
% Full 2CXM Model, no Delay       (4 parameters, Flow, Ktrans, Vp, Ve)
Init_Guess_Larsson_no_Delay               = double( [est_Vb_Patlak_noise_no_Delay est_E_Patlak_noise_no_Delay init_Ve_guess] );
est_F_no_Delay                            = double(Flow_no_Delay);

% Highly Perfused, F->Inf + Delay (4 parameters, Ktrans, Vp, Ve, Delay) - Same as ETM
Init_Guess_Larsson_with_Delay_High_F      = double( [est_Vb_Patlak_noise_with_Delay est_Ktrans_Patlak_noise_with_Delay init_Ve_guess] );
est_F_with_Delay_High_F                   = double(exp(100)); % Infinity
% Highly Perfused, F->Inf         (3 parameters, Ktrans, Vp, Ve) - Same as ETM
Init_Guess_Larsson_no_Delay_High_F        = double( [est_Vb_Patlak_noise_no_Delay est_Ktrans_Patlak_noise_no_Delay init_Ve_guess] );
est_F_no_Delay_High_F                     = double(exp(100)); % Infinity

% Highly Perfused, Ve->0 + Delay (4 parameters, Ktrans, Vp, Ve, Delay) - Same as ETM
Init_Guess_Larsson_with_Delay_no_Ve       = double( [est_Vb_Patlak_noise_with_Delay est_E_Patlak_noise_with_Delay ] );
est_F_with_Delay_no_Ve                    = double(Flow_with_Delay);
% Highly Perfused, Ve->0         (3 parameters, Ktrans, Vp, Ve) - Same as ETM
Init_Guess_Larsson_no_Delay_no_Ve         = double( [est_Vb_Patlak_noise_no_Delay est_E_Patlak_noise_no_Delay ] );
est_F_no_Delay_no_Ve                      = double(Flow_no_Delay);

% No Indicator Exchange/Highly Vascularized  : Ve->0, E->0 +Delay (3 parameters, Flow, Vp, Delay)
Init_Guess_Larsson_with_Delay_no_E        = double( est_Vb_Patlak_noise_with_Delay );
est_F_with_Delay_no_E                     = double(Flow_with_Delay);
% No Indicator Exchange/Highly Vascularized  : Ve->0, E->0 (2 parameters, Flow, Vp)
Init_Guess_Larsson_no_Delay_no_E          = double( est_Vb_Patlak_noise_no_Delay );
est_F_no_Delay_no_E                       = double(Flow_no_Delay);

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
        
    Larsson_function_with_Delay_no_Ve       = @(x,t) Adjusted_Larsson_Filter_no_Ve       ( t, est_F_with_Delay      , x(1), x(2)      );
    Larsson_function_no_Delay_no_Ve         = @(x,t) Adjusted_Larsson_Filter_no_Ve       ( t, est_F_no_Delay        , x(1), x(2)      );
    
    Larsson_function_with_Delay_High_F      = @(x,t) Adjusted_Larsson_Filter_High_F     ( t,                         x(1), x(2), x(3));
    Larsson_function_no_Delay_High_F        = @(x,t) Adjusted_Larsson_Filter_High_F     ( t,                         x(1), x(2), x(3));
    
    Larsson_function_with_Delay_no_E        = @(x,t) Adjusted_Larsson_Filter_no_E       ( t, est_F_with_Delay_no_E , x(1)            );
    Larsson_function_no_Delay_no_E          = @(x,t) Adjusted_Larsson_Filter_no_E       ( t, est_F_no_Delay_no_E   , x(1)            );
    
    Larsson_function_with_Delay_no_E_High_F = @(x,t) Adjusted_Larsson_Filter_no_E_High_F( t,                         x(1)            );
    Larsson_function_no_Delay_no_E_High_F   = @(x,t) Adjusted_Larsson_Filter_no_E_High_F( t,                         x(1)            );
else
    Larsson_function_with_Delay             = @(x,t) Larsson_Filter                     ( t, est_F_with_Delay             , x(1), x(2), x(3), Hct);
    Larsson_function_no_Delay               = @(x,t) Larsson_Filter                     ( t, est_F_no_Delay               , x(1), x(2), x(3), Hct);
    
    Larsson_function_with_Delay_High_F      = @(x,t) Larsson_Filter_High_F              ( t,                                x(1), x(2), x(3), Hct);
    Larsson_function_no_Delay_High_F        = @(x,t) Larsson_Filter_High_F              ( t,                                x(1), x(2), x(3), Hct);
    
    Larsson_function_with_Delay_no_E        = @(x,t) Larsson_Filter_no_E                ( t, est_F_with_Delay_no_E        , x(1)            , Hct);
    Larsson_function_no_Delay_no_E          = @(x,t) Larsson_Filter_no_E                ( t, est_F_no_Delay_no_E          , x(1)            , Hct);
    
    Larsson_function_with_Delay_no_E_High_F = @(x,t) Larsson_Filter_no_E_High_F         ( t,                                x(1)            , Hct);
    Larsson_function_no_Delay_no_E_High_F   = @(x,t) Larsson_Filter_no_E_High_F         ( t,                                x(1)            , Hct);
end


% Calculate parameters only in case F is not 0 (or else there is a problem)
if (est_F_no_Delay~=0)
    
    if Correct_estimation_due_to_delay
        % Fitting with delay
        %[est_params_Larsson_with_Delay_noise,residue_norm_Larsson_with_Delay_noise,residual_Larsson_with_Delay_noise,exitflag_Larsson_with_Delay_noise,algo_info_Larsson_with_Delay_noise] = ...
        [est_params_Larsson_with_Delay_noise,~,~,~,~] = ...
            lsqcurvefit(Larsson_function_with_Delay,Init_Guess_Larsson_with_Delay,time_vec_minutes_T,Est_ht_with_Delay_T/Flow_with_Delay,...
            LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
        
        % Fitting with no delay
        if ~Ignore_Delay_Model_Selection 
            %[est_params_Larsson_no_Delay_noise,residue_norm_Larsson_no_Delay_noise,residual_Larsson_no_Delay_noise,exitflag_Larsson_no_Delay_noise,algo_info_Larsson_no_Delay_noise] = ...
            [est_params_Larsson_no_Delay_noise,~,~,~,~] = ...
                lsqcurvefit(Larsson_function_no_Delay,Init_Guess_Larsson_no_Delay,time_vec_minutes_T,Est_ht_no_Delay_T/Flow_no_Delay,...
                LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
        end        
    else
        % With and without delay is the same
        
        % Fitting with delay
        %[est_params_Larsson_with_Delay_noise,residue_norm_Larsson_with_Delay_noise,residual_Larsson_with_Delay_noise,exitflag_Larsson_with_Delay_noise,algo_info_Larsson_with_Delay_noise] = ...
        [est_params_Larsson_with_Delay_noise,~,~,~,~] = ...
            lsqcurvefit(Larsson_function_with_Delay,Init_Guess_Larsson_with_Delay,time_vec_minutes_T,Est_ht_with_Delay_T/Flow_with_Delay,...
            LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
        
        % Fitting with no delay
        est_params_Larsson_no_Delay_noise   = est_params_Larsson_with_Delay_noise;
    end
    
    % Assigning two compartment parameters estimation
    Vb_with_Delay         = est_params_Larsson_with_Delay_noise(1);
    Vb_no_Delay           = est_params_Larsson_no_Delay_noise(1);
    E_with_Delay          = est_params_Larsson_with_Delay_noise(2);
    E_no_Delay            = est_params_Larsson_no_Delay_noise(2);
    Ktrans_with_Delay     = est_params_Larsson_with_Delay_noise(2)* Flow_with_Delay;
    Ktrans_no_Delay       = est_params_Larsson_no_Delay_noise(2)* Flow_no_Delay;
    Ve_with_Delay         = est_params_Larsson_with_Delay_noise(3);
    Ve_no_Delay           = est_params_Larsson_no_Delay_noise(3);
    
    % Estimate MTT
    est_IRF_with_Delay           = Est_ht_with_Delay_T / Flow_with_Delay;
    est_IRF_no_Delay             = Est_ht_no_Delay_T  / Flow_no_Delay;
    
    est_MTT_noise_with_Delay     = cumtrapz(time_vec_minutes,est_IRF_with_Delay);
    MTT_with_Delay               = est_MTT_noise_with_Delay(end);
    
    est_MTT_noise_no_Delay       = cumtrapz(time_vec_minutes,est_IRF_no_Delay);
    MTT_no_Delay                 = est_MTT_noise_no_Delay(end);
    
    if (Adjusted_Larsson_Model)
        fitted_larsson_with_Delay = Flow_with_Delay*Adjusted_Larsson_Filter( time_vec_minutes, Flow_with_Delay, Vb_with_Delay, E_with_Delay, Ve_with_Delay);
        fitted_larsson_no_Delay   = Flow_no_Delay  *Adjusted_Larsson_Filter( time_vec_minutes, Flow_no_Delay, Vb_no_Delay, E_no_Delay, Ve_no_Delay);
    else
        fitted_larsson_with_Delay = Flow_with_Delay*Larsson_Filter( time_vec_minutes, Flow_with_Delay, Vb_with_Delay, E_with_Delay, Ve_with_Delay, Hct);
        fitted_larsson_no_Delay   = Flow_no_Delay  *Larsson_Filter( time_vec_minutes, Flow_no_Delay, Vb_no_Delay, E_no_Delay, Ve_no_Delay, Hct);
    end
    
    
    %% Model Selection
    if Use_Model_Selection
        
        % ------------------------------------- Fit with no Ve + Delay ---------------------------------------------------------------
        %To_Fit_High_F_with_Delay = Est_ht_with_Delay_T/Flow_with_Delay;
        To_Fit_no_Ve_with_Delay = Est_ht_with_Delay_T;
        %[est_params_Larsson_noise_with_Delay_High_F,residue_norm_Larsson_noise_with_Delay_High_F,residual_Larsson_noise_with_Delay_High_F,exitflag_Larsson_noise_with_Delay_High_F,algo_info_Larsson_noise_with_Delay_High_F] = ...
        [est_params_Larsson_noise_with_Delay_no_Ve,~,~,~,~] = ...
            lsqcurvefit(Larsson_function_with_Delay_no_Ve,Init_Guess_Larsson_with_Delay_no_Ve,time_vec_minutes_T,To_Fit_no_Ve_with_Delay,...
            LowerBound_Larsson_no_Ve,UpperBound_Larsson_no_Ve,algorithm_options);
        
        % Assigning two compartment parameters estimation
        Vb_with_Delay_no_Ve       = est_params_Larsson_noise_with_Delay_no_Ve(1);
        E_with_Delay_no_Ve        = est_params_Larsson_noise_with_Delay_no_Ve(2);
        
        % ------------------------------------- Fit with no Ve ---------------------------------------------------------------
        if ~Ignore_Delay_Model_Selection
            
            %To_Fit_High_F_no_Delay = Est_ht_with_Delay_T/Flow_no_Delay;
            To_Fit_no_Ve_no_Delay = Est_ht_no_Delay_T;
            %[est_params_Larsson_noise_no_Delay_High_F,residue_norm_Larsson_noise_no_Delay_High_F,residual_Larsson_noise_no_Delay_High_F,exitflag_Larsson_noise_no_Delay_High_F,algo_info_Larsson_noise_no_Delay_High_F] = ...
            [est_params_Larsson_noise_no_Delay_no_Ve,~,~,~,~] = ...
                lsqcurvefit(Larsson_function_no_Delay_no_Ve,Init_Guess_Larsson_no_Delay_no_Ve,time_vec_minutes_T,To_Fit_no_Ve_no_Delay,...
                LowerBound_Larsson_no_Ve,UpperBound_Larsson_no_Ve,algorithm_options);
            
            % Assigning two compartment parameters estimation
            Vb_no_Delay_no_Ve       = est_params_Larsson_noise_no_Delay_no_Ve(1);
            E_no_Delay_no_Ve        = est_params_Larsson_noise_no_Delay_no_Ve(2);
            
        end
        
        % ------------------------------------- Fit with high f + Delay ---------------------------------------------------------------
        %To_Fit_High_F_with_Delay = Est_ht_with_Delay_T/Flow_with_Delay;
        To_Fit_High_F_with_Delay = Est_ht_with_Delay_T;
        %[est_params_Larsson_noise_with_Delay_High_F,residue_norm_Larsson_noise_with_Delay_High_F,residual_Larsson_noise_with_Delay_High_F,exitflag_Larsson_noise_with_Delay_High_F,algo_info_Larsson_noise_with_Delay_High_F] = ...
        [est_params_Larsson_noise_with_Delay_High_F,~,~,~,~] = ...
            lsqcurvefit(Larsson_function_with_Delay_High_F,Init_Guess_Larsson_with_Delay_High_F,time_vec_minutes_T,To_Fit_High_F_with_Delay,...
            LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
        
        % Assigning two compartment parameters estimation
        Vb_with_Delay_High_F       = est_params_Larsson_noise_with_Delay_High_F(1);
        Ktrans_with_Delay_High_F   = est_params_Larsson_noise_with_Delay_High_F(2);
        Ve_with_Delay_High_F       = est_params_Larsson_noise_with_Delay_High_F(3);
        
        % ------------------------------------- Fit with high f ---------------------------------------------------------------
        if ~Ignore_Delay_Model_Selection
            
            %To_Fit_High_F_no_Delay = Est_ht_with_Delay_T/Flow_no_Delay;
            To_Fit_High_F_no_Delay = Est_ht_no_Delay_T;
            %[est_params_Larsson_noise_no_Delay_High_F,residue_norm_Larsson_noise_no_Delay_High_F,residual_Larsson_noise_no_Delay_High_F,exitflag_Larsson_noise_no_Delay_High_F,algo_info_Larsson_noise_no_Delay_High_F] = ...
            [est_params_Larsson_noise_no_Delay_High_F,~,~,~,~] = ...
                lsqcurvefit(Larsson_function_no_Delay_High_F,Init_Guess_Larsson_no_Delay_High_F,time_vec_minutes_T,To_Fit_High_F_no_Delay,...
                LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
            
            % Assigning two compartment parameters estimation
            Vb_no_Delay_High_F       = est_params_Larsson_noise_no_Delay_High_F(1);
            Ktrans_no_Delay_High_F   = est_params_Larsson_noise_no_Delay_High_F(2);
            Ve_no_Delay_High_F       = est_params_Larsson_noise_no_Delay_High_F(3);
            
        end
        % ------------------------------------- Fit with delay and no permeability ---------------------------------------------------------------
        To_Fit_no_E_with_Delay = Est_ht_with_Delay_T/Flow_with_Delay;
        
        %[est_params_Larsson_noise_with_Delay_no_E,residue_norm_Larsson_noise_with_Delay_no_E,residual_Larsson_noise_with_Delay_no_E,exitflag_Larsson_noise_with_Delay_no_E,algo_info_Larsson_noise_with_Delay_no_E] = ...
        [est_params_Larsson_noise_with_Delay_no_E,~,~,~,~] = ...
            lsqcurvefit(Larsson_function_with_Delay_no_E,Init_Guess_Larsson_with_Delay_no_E,time_vec_minutes_T,To_Fit_no_E_with_Delay,...
            LowerBound_Larsson_no_E,UpperBound_Larsson_no_E,algorithm_options);
        
        % Assigning two compartment parameters estimation
        Vb_with_Delay_no_E         = est_params_Larsson_noise_with_Delay_no_E(1);
        
        % ------------------------------------- Fit with no delay and no permeability ---------------------------------------------------------------
        if ~Ignore_Delay_Model_Selection
            
            To_Fit_no_E_no_Delay = Est_ht_no_Delay_T/Flow_no_Delay;
            
            %[est_params_Larsson_noise_no_Delay_no_E,residue_norm_Larsson_noise_no_Delay_no_E,residual_Larsson_noise_no_Delay_no_E,exitflag_Larsson_noise_no_Delay_no_E,algo_info_Larsson_noise_no_Delay_no_E] = ...
            [est_params_Larsson_noise_no_Delay_no_E,~,~,~,~] = ...
                lsqcurvefit(Larsson_function_no_Delay_no_E,Init_Guess_Larsson_no_Delay_no_E,time_vec_minutes_T,To_Fit_no_E_no_Delay,...
                LowerBound_Larsson_no_E,UpperBound_Larsson_no_E,algorithm_options);
            
            % Assigning two compartment parameters estimation
            Vb_no_Delay_no_E         = est_params_Larsson_noise_no_Delay_no_E(1);
            
        end
        
        % ------------------------------------- Fit with high f and no permeability + Delay ---------------------------------------------------------------
        %To_Fit_no_E_High_F_with_Delay = Est_ht_with_Delay_T/Flow_with_Delay;
        To_Fit_no_E_High_F_with_Delay = Est_ht_with_Delay_T;
        %[est_params_Larsson_noise_with_Delay_no_E_High_F,residue_norm_Larsson_noise_with_Delay_no_E_High_F,residual_Larsson_noise_with_Delay_no_E_High_F,exitflag_Larsson_noise_with_Delay_no_E_High_F,algo_info_Larsson_noise_with_Delay_no_E_High_F] = ...
        [est_params_Larsson_noise_with_Delay_no_E_High_F,~,~,~,~] = ...
            lsqcurvefit(Larsson_function_with_Delay_no_E_High_F,Init_Guess_Larsson_with_Delay_no_E_High_F,time_vec_minutes_T,To_Fit_no_E_High_F_with_Delay,...
            LowerBound_Larsson_no_E,UpperBound_Larsson_no_E,algorithm_options);
        
        % Assigning two compartment parameters estimation
        Vb_with_Delay_no_E_High_F         = est_params_Larsson_noise_with_Delay_no_E_High_F(1);
        
        % ------------------------------------- Fit with high f and no permeability ---------------------------------------------------------------
        %To_Fit_no_E_High_F_no_Delay = Est_ht_with_Delay_T/Flow_no_Delay;
        To_Fit_no_E_High_F_no_Delay = Est_ht_no_Delay_T;
        if ~Ignore_Delay_Model_Selection
            %[est_params_Larsson_noise_no_Delay_no_E_High_F,residue_norm_Larsson_noise_no_Delay_no_E_High_F,residual_Larsson_noise_no_Delay_no_E_High_F,exitflag_Larsson_noise_no_Delay_no_E_High_F,algo_info_Larsson_noise_no_Delay_no_E_High_F] = ...
            [est_params_Larsson_noise_no_Delay_no_E_High_F,~,~,~,~] = ...
                lsqcurvefit(Larsson_function_no_Delay_no_E_High_F,Init_Guess_Larsson_no_Delay_no_E_High_F,time_vec_minutes_T,To_Fit_no_E_High_F_no_Delay,...
                LowerBound_Larsson_no_E,UpperBound_Larsson_no_E,algorithm_options);
            
            % Assigning two compartment parameters estimation
            Vb_no_Delay_no_E_High_F         = est_params_Larsson_noise_no_Delay_no_E_High_F(1);
        end
                
        
        if (Adjusted_Larsson_Model)
            
            % Effectively ETM model - 3 parameters
            fitted_larsson_with_Delay_High_F      = Adjusted_Larsson_Filter_High_F( time_vec_minutes, Vb_with_Delay_High_F , Ktrans_with_Delay_High_F , Ve_with_Delay_High_F);
            if ~Ignore_Delay_Model_Selection
                fitted_larsson_no_Delay_High_F    = Adjusted_Larsson_Filter_High_F( time_vec_minutes, Vb_no_Delay_High_F   , Ktrans_no_Delay_High_F   , Ve_no_Delay_High_F);
            end
            
            % Uptake model -> should fit healthy brain tissues - 3 parameters
            fitted_larsson_with_Delay_no_Ve      = Adjusted_Larsson_Filter_no_Ve( time_vec_minutes, Flow_with_Delay, Vb_with_Delay_no_Ve , E_with_Delay_no_Ve);
            if ~Ignore_Delay_Model_Selection
                fitted_larsson_no_Delay_no_Ve    = Adjusted_Larsson_Filter_no_Ve( time_vec_minutes, Flow_no_Delay, Vb_no_Delay_no_Ve  , E_no_Delay_no_Ve);
            end
            
            % Uptake model - no permeability -> should fit WM
            fitted_larsson_with_Delay_no_E            = Flow_with_Delay * Adjusted_Larsson_Filter_no_E( time_vec_minutes, Flow_with_Delay , Vb_with_Delay_no_E);
            if ~Ignore_Delay_Model_Selection
                fitted_larsson_no_Delay_no_E          = Flow_no_Delay   * Adjusted_Larsson_Filter_no_E( time_vec_minutes, Flow_no_Delay   , Vb_no_Delay_no_E);
            end
            
            % Impulse response at height Vp -> should fit GM
            fitted_larsson_with_Delay_no_E_High_F = Adjusted_Larsson_Filter_no_E_High_F ( time_vec_minutes,                          Vb_with_Delay_no_E_High_F );
            if ~Ignore_Delay_Model_Selection
                fitted_larsson_no_Delay_no_E_High_F   = Adjusted_Larsson_Filter_no_E_High_F ( time_vec_minutes,                          Vb_no_Delay_no_E_High_F   );
            end
        else
            
            fitted_larsson_with_Delay_High_F      = Larsson_function_with_Delay_High_F( time_vec_minutes, est_F_with_Delay_High_F      , Vb_with_Delay_High_F      , 0, Ve_with_Delay_High_F, Hct);
            fitted_larsson_no_Delay_High_F        = Larsson_function_no_Delay_High_F( time_vec_minutes, est_F_no_Delay_High_F        , Vb_no_Delay_High_F        , 0, Ve_no_Delay_High_F  , Hct);
            
            fitted_larsson_with_Delay_no_Ve       = Larsson_function_with_Delay_no_Ve ( time_vec_minutes, est_F_with_Delay_no_Ve       , Vb_with_Delay_no_Ve       , E_with_Delay_no_Ve, Hct);
            fitted_larsson_no_Delay_no_Ve         = Larsson_function_no_Delay_no_Ve ( time_vec_minutes, est_F_no_Delay_no_Ve         , Vb_no_Delay_no_Ve         , E_no_Delay_no_Ve  , Hct);
            
            fitted_larsson_with_Delay_no_E        = Flow_with_Delay * Larsson_function_with_Delay_no_E( time_vec_minutes, Flow_with_Delay       , Vb_with_Delay_no_E        , 0, 0                          , Hct);
            fitted_larsson_no_Delay_no_E          = Flow_no_Delay   * Larsson_function_no_Delay_no_E( time_vec_minutes, Flow_no_Delay         , Vb_no_Delay_no_E          , 0, 0                          , Hct);
            fitted_larsson_with_Delay_no_E_High_F = Larsson_function_with_Delay_no_E_High_F( time_vec_minutes, est_F_with_Delay_no_E_High_F , Vb_with_Delay_no_E_High_F , 0, 0                          , Hct);
            fitted_larsson_no_Delay_no_E_High_F   = Larsson_function_no_Delay_no_E_High_F( time_vec_minutes, est_F_no_Delay_no_E_High_F   , Vb_no_Delay_no_E_High_F   , 0, 0                          , Hct);
        end
        
    end % if use_model_selection
    
    % Assigning Patlak parameters estimation
    Ktrans_with_Delay_Patlak    = est_Ktrans_Patlak_noise_with_Delay;
    Ktrans_no_Delay_Patlak      = est_Ktrans_Patlak_noise_no_Delay;
    Vb_with_Delay_Patlak        = est_Vb_Patlak_noise_with_Delay;
    Vb_no_Delay_Patlak          = est_Vb_Patlak_noise_no_Delay;
    MTT_with_Delay_Patlak       = est_MTT_Patlak_noise_with_Delay;
    MTT_no_Delay_Patlak         = est_MTT_Patlak_noise_no_Delay;
    
    % Delay of the AIF will be calculated according to the place of the maximum value of F*IRF
    max_index_with_Delay               = find( Est_ht_with_Delay_T == Flow_with_Delay, 1);
    max_index_no_Delay                 = find( Est_ht_no_Delay_T   == Flow_no_Delay, 1);
    
    % Translate to minutes
    Delay_sec_by_Max_Val_with_Delay = (max_index_with_Delay - 1) * min_interval * 60;
    Delay_sec_by_Max_Val_no_Delay   = (max_index_no_Delay   - 1) * min_interval * 60;
    
end

%% Return Parameters
% SINGLE GAUSSIAN: 
% -----------------------------
% gaussian_param
% fitted_gaussian

% DOUBLE GAUSSIAN:
% -----------------------------
% double_gaussian_param
% fitted_double_gaussian

% LARSSON
% -----------------------------
% Flow_with_Delay
% Flow_no_Delay
% Vb_with_Delay
% Vb_no_Delay
% E_with_Delay
% E_no_Delay
% Ktrans_with_Delay
% Ktrans_no_Delay
% Ve_with_Delay
% Ve_no_Delay
% MTT_with_Delay
% MTT_no_Delay
% fitted_larsson_with_Delay
% fitted_larsson_no_Delay

% Model Selection - Fit with high f + Delay
% -----------------------------
% Vb_with_Delay_High_F
% Ktrans_with_Delay_High_F
% Ve_with_Delay_High_F
% fitted_larsson_with_Delay_High_F

% Model Selection -  Fit with high f
% -----------------------------
% Vb_no_Delay_High_F
% Ktrans_no_Delay_High_F
% Ve_no_Delay_High_F
% fitted_larsson_no_Delay_High_F

% Model Selection - Fit with no Ve + Delay (Uptake)
% -----------------------------
% Vb_with_Delay_no_Ve
% E_with_Delay_no_Ve
% fitted_larsson_with_Delay_no_Ve

% Model Selection -  Fit with no Ve       (Uptake)
% -----------------------------
% Vb_no_Delay_no_Ve
% E_no_Delay_no_Ve
% fitted_larsson_no_Delay_no_Ve

% Model Selection - Fit with delay and no permeability
% -----------------------------
% Vb_with_Delay_no_E
% fitted_larsson_with_Delay_no_E

% Model Selection - Fit with no delay and no permeability
% -----------------------------
% Vb_no_Delay_no_E
% fitted_larsson_no_Delay_no_E

% Model Selection - Fit with high f and no permeability + Delay
% -----------------------------
% Vb_with_Delay_no_E_High_F
% fitted_larsson_with_Delay_no_E_High_F

% Model Selection - Fit with high f and no permeability
% -----------------------------
% Vb_no_Delay_no_E_High_F
% fitted_larsson_no_Delay_no_E_High_F

% Patlak
% -----------------------------
% Ktrans_with_Delay_Patlak
% Ktrans_no_Delay_Patlak
% Vb_with_Delay_Patlak
% Vb_no_Delay_Patlak
% MTT_with_Delay_Patlak
% MTT_no_Delay_Patlak
% Delay_sec_by_Max_Val_with_Delay
% Delay_sec_by_Max_Val_no_Delay

end

