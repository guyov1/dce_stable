function [ Sim_Struct ] = Simulation_Set_Params( Sim_Struct, Verbosity )

if strcmp(Verbosity,'Full')
    display('-I- Setting simulation parameters...');
end

% Initiate index for future figures
Sim_Struct.idx_fig = 1;

%% AIF parameters - Parker's AIF

% Factor for AIF to have realistic values like we see in MRI images
Sim_Struct.r_factor = 0.0116;

% AIF average population parameters
Sim_Struct.A1    = 0.809;
Sim_Struct.A2    = 0.330;

% Changed to shift the boluses to the right
Sim_Struct.T1    = 0.17046;
Sim_Struct.T2    = 0.365;
%T1    = 0.27046;
%T2    = 0.465;

Sim_Struct.sig1  = 0.0563;
Sim_Struct.sig2  = 0.132;
Sim_Struct.alpha = 1.050;
Sim_Struct.beta  = 0.1685;
Sim_Struct.s     = 38.078;
Sim_Struct.tau   = 0.483;

%%

% Set number of iterations for simulation
Sim_Struct.num_iterations = 15;
% Do each iteration a few time and average results for better statistic information
Sim_Struct.num_averages   = 5;

% Gaussian filter parameters
Sim_Struct.sigma                        = 2/60;   % In minutes
Sim_Struct.sigma_vec                    = linspace(0.5/60,50/60,Sim_Struct.num_iterations); % When iterating
Sim_Struct.t_d                          = 30/60;  % Time shift in minutes
Sim_Struct.t_d_vec                      = linspace(5/60,50/60,Sim_Struct.num_iterations); % When iterating;
Sim_Struct.amplitude                    = 1;      % Amplitude of Gaussian
Sim_Struct.amplitude_vec                = linspace(0.1,5,Sim_Struct.num_iterations); % When iterating;;

% Number of FFT points
Sim_Struct.num_FFT_points = 2000;

% Total simulation time
Sim_Struct.total_sim_time_min  = 4; %[min]

% Non-linear square fitting parameters
Sim_Struct.FMS_TolFun          = 1e-11;
Sim_Struct.FMS_MaxFunEvals     = 10000;
Sim_Struct.FMS_MaxIter         = 10000;
%FMS_Algorithm      = 'levenberg-marquardt'; % Does not supper lower-upper bounds
Sim_Struct.FMS_Algorithm       = 'trust-region-reflective';
Sim_Struct.algorithm_options   = optimset('TolFun',Sim_Struct.FMS_TolFun,'MaxFunEvals',Sim_Struct.FMS_MaxFunEvals,'MaxIter',Sim_Struct.FMS_MaxIter,...
    'Display','off','Algorithm',Sim_Struct.FMS_Algorithm);

% Initial Guess for non-linear curve fitting for Gaussian (td=1sec, var=1[sec], amplitude = 1)
Sim_Struct.Init_Guess_Gaussian = [1/60 1/60 0.5];

% Set lower and upper bounds for parameters for non-linear fitting
% Gaussian

% Time delay -> 0 to 60 seconds
% Sigma      -> 0 to 30 seconds
% Amplitude  -> 0 to 3
Sim_Struct.LowerBound_Gauss  = [0     0     0];
Sim_Struct.UpperBound_Gauss  = [60/60 30/60 3];

% Larsson parameters boundaries for non-linear curve fitting
% Vb -> 0 to 100 [ml/100g]
% E  -> 0 to 1
% Ve -> 0 to 100 [ml/100g]
Sim_Struct.LowerBound_Larsson  = [0   0     0];
Sim_Struct.UpperBound_Larsson  = [100 1     100];

% Additional AIF parameters
Sim_Struct.additional_AIF_delay_sec     = 0.0; % Delay added to AIF before filtering
Sim_Struct.additional_AIF_delay_sec_vec = linspace(-2,3,Sim_Struct.num_iterations); % When iterating;;

% Larsson filter parameters
Sim_Struct.F_single   = 60;                                         % When a single iteration
Sim_Struct.F_vec      = linspace(10,140,Sim_Struct.num_iterations); % When iterating
Sim_Struct.Vb_single  = 12;
Sim_Struct.Vb_vec     = linspace(3,20,Sim_Struct.num_iterations);
Sim_Struct.E_single   = 0.1;
Sim_Struct.E_vec      = linspace(0.01,0.99,Sim_Struct.num_iterations);
Sim_Struct.Ve_single  = 0.05; % 0.1
Sim_Struct.Ve_vec     = linspace(0,1,Sim_Struct.num_iterations);
Sim_Struct.Hct_single = 0.38;

% Determines SNR ( noise_var = mean(signal)/SNR_base )
Sim_Struct.SNR_single = 15;
Sim_Struct.SNR_vec    = linspace(20,1,Sim_Struct.num_iterations);

% Time interval between samples
Sim_Struct.sec_interval     = 2;                                          %[sec]
Sim_Struct.sec_vec          = linspace(0.5,15,Sim_Struct.num_iterations); % For iterations
% Round to the nearest 0.5
Sim_Struct.sec_vec          = round(Sim_Struct.sec_vec*2)/2;
Sim_Struct.Fs               = 1 / Sim_Struct.sec_interval; %[Hz] - Sampling rate
Sim_Struct.min_interval     = Sim_Struct.sec_interval/60;  %[min]
Sim_Struct.num_time_stamps  = round(Sim_Struct.total_sim_time_min / Sim_Struct.min_interval);
% Time vector for AIF and Ct(t)
Sim_Struct.time_vec_minutes = (0:Sim_Struct.num_time_stamps-1).* Sim_Struct.min_interval;
% High resolution before downsampling to min_interval
Sim_Struct.High_res_min     = 0.01/60;

% Determine according to what parameter to check simulations
Sim_Struct.iterate_SNR                   = 0;
Sim_Struct.iterate_sec_interval          = 0;
Sim_Struct.iterate_gaussian_sigma        = 0;
Sim_Struct.iterate_gaussian_time_delay   = 0;
Sim_Struct.iterate_gaussian_amplitude    = 0;
Sim_Struct.iterate_F_larsson             = 1;
Sim_Struct.iterate_Vb_larsson            = 0;
Sim_Struct.iterate_E_larsson             = 0;
Sim_Struct.iterate_Ve_larsson            = 0;
Sim_Struct.iterate_AIF_delay             = 0;

% If not iterating, update back to 1 iteration
if ~(Sim_Struct.iterate_SNR || Sim_Struct.iterate_sec_interval || Sim_Struct.iterate_gaussian_sigma || Sim_Struct.iterate_gaussian_time_delay || ...
        Sim_Struct.iterate_gaussian_amplitude || Sim_Struct.iterate_F_larsson || Sim_Struct.iterate_Vb_larsson || Sim_Struct.iterate_E_larsson || Sim_Struct.iterate_Ve_larsson || Sim_Struct.iterate_AIF_delay)
    
    Sim_Struct.num_iterations = 1;
    Sim_Struct.num_averages   = 1;
    
end

% Plot graphs when not iterating
% When iterating, plot the error as a function of parameters
if (Sim_Struct.num_iterations == 1 && Sim_Struct.num_averages == 1)
    Sim_Struct.plot_flag = true;
    Sim_Struct.plot_error_results_flag = false;
else
    Sim_Struct.plot_flag = false;
    Sim_Struct.plot_error_results_flag = true;
end

% Additional parameters for simulation iterations
Sim_Struct.Ignore_Gaussian_Calculation   = 1;      % Ignores all calculations relating to gaussian function
Sim_Struct.Iterate_Murase_Tofts          = 0;      % Check Murase calculation for Tofts parameters
Sim_Struct.Iterate_Murase_Tofts_num_iter = 100000; % number of iterations for Murase Tofts check
Sim_Struct.One_Iteration_Murase_Tofts    = 0;

% Check Sourbron 4 parameters estimation
Sim_Struct.Check_Sourbron_Estimate       = 1;
Sim_Struct.Random_init_F_guess           = false; % If false, take larsson F estimation

% Drive differential equation
Sim_Struct.Drive_Diff_Eq                 = 0;

% If checking Murase vs. Tofts, avoid averages iterations
if (Sim_Struct.Iterate_Murase_Tofts)
    Sim_Struct.num_averages = 1;
end

% In case we want one iteration for graph display, change # to 1
if (Sim_Struct.One_Iteration_Murase_Tofts)
    Sim_Struct.num_iterations = 1;
    Sim_Struct.num_averages   = 1;
end

% Set if derivative matrix will be divided by time
Sim_Struct.Derivative_Time_Devision = false;

% Choose which Patlak estimation to take
% Possible - 1. "Specified Points" 2. "All Points" 3. "Weighted Points"
Sim_Struct.Patlak_Est = 'Specified Points';

% Choose which filter estimation to use when calculating parameters (non-linear method)
% Possible: 'Wiener', 'Ridge', 'Spline', 'Spline_1st', 'Spline_2nd', 'PCA'
Sim_Struct.Filter_Est_Chosen = 'PCA';

% Apply cyclic convolution to compensate for AIF delay
Sim_Struct.Use_Cyclic_Conv_4_ht_est               = false;     % Use cyclic de-convolution to correct for delay
Sim_Struct.Use_Upsampling_and_Cyclic              = false;     % Use cyclic de-convolution to correct for delay + upsampling
Sim_Struct.Use_Upsampling_Delay_Comp              = false;     % Upsample Ct(t) and AIF(t) to try and predict time shift in AIF
Sim_Struct.Upsampling_resolution                  = 0.1 / 60;  % Set the upsampling target
Sim_Struct.Correct_estimation_due_to_delay        = false;      % Try to correct for delay
Sim_Struct.RMS_Smooth                             = true;      % When calculating RMS, smooth CTC first
Sim_Struct.RMS_Smooth_Around_Bolus                = false;     % Calculate RMS around bolus only (to avoid noise aggregation afterwards)
Sim_Struct.Simple_AIF_Delay_Correct               = false;     % Correct AIF by max point shift
Sim_Struct.Max_Time_Delay                         = 3;         % Set the maximal possible time delay in seconds
Sim_Struct.Diff_From_Bolus                        = 10;        % The difference in seconds from the bolus to look on
Sim_Struct.BiExp2CTC_RMS_Ratio                    = 0;         % Sets the ratio between BiExp fit and CTC fit when estimating time delay

% Polynomial degree for basis splines
Sim_Struct.poly_deg      = 4;

% Choose knots for splines (currently takes every 1 out of 2 points)
Sim_Struct.knot_interval = 5;
Sim_Struct.knots         = Sim_Struct.time_vec_minutes(1:Sim_Struct.knot_interval:end);

% Initiate with zeros the estimated delay and shifted AIF which are later calculated
Sim_Struct.est_delay_by_AIF_correct         = zeros(Sim_Struct.num_time_stamps,1);
Sim_Struct.Sim_AIF_with_noise_Regul_shifted = zeros(Sim_Struct.num_time_stamps,1);

% Lambda for regularization
% 1. Ridge regression  2. Spline-no deriv
% 3. Spline-1st deriv  4. Spline-2nd deriv
%lambda_vec              = [9 132 0.0391 0.0141];
Sim_Struct.lambda_vec_gauss        = [34.7 11 0.1 16.8];
%lambda_vec              = [9 132 0.0391 0.4];
Sim_Struct.lambda_vec_larss        = [4.7 0.9 0.1 0.4];
%Sim_Struct.lambda_vec_larss        = [9 132 0.0391 1.9];
%Sim_Struct.lambda_vec_larss        = [9 132 0.0391 0.012];
Sim_Struct.lambda_vec_larss        = [4.7 0.9 0.1 0.2];
%lambda_vec              = [9 132 0.0391 1.9];
%lambda_vec              = [9 132 0.0391 0.0];
%lambda_vec              = [9 132 0.0391 0.012]; 74.8
Sim_Struct.normalize               = 1;  % Normalize flag for ridge() function

% Choose whether to plot L curve
Sim_Struct.plot_L_Curve = 0;

% Initiate results matrix which will hold:
% 1.     SNR ratio
% 2.     Seconds interval of simulation
% 3-6.   Gaussian:        original sigma,      estimated sigma, error percent and standard deviation
% 7-10.  Gaussian:        original delay time, estimated delay time, error percent and standard deviation
% 11-14. Gaussian:        original amplitude,  estimated amplitude, error percent and standard deviation
% 15-18. Larsson filter:  original Flow,       estimated Flow, error percent and standard deviation
% 19-22. AIF:             original delay,      estimated delay, error percent and standard deviation
% 23-26. Larsson filter:  original Ki,         estimated Ki using Patlak, error percent and standard deviation
% 27-30. Larsson filter:  original Ki,         estimated Ki using 2CXM, error percent and standard deviation
% 31-34. Larsson filter:  original PS,         estimated PS, error percent and standard deviation
% 35-38. Larsson filter:  original Vb,         estimated Vb using Patlak, error percent and standard deviation
% 39-42. Larsson filter:  original Vb,         estimated Vb using 2CXM, error percent and standard deviation
% 43-46. Larsson filter:  original Vd,         estimated Vd, error percent and standard deviation
% 47-50. Larsson filter:  original Vd,         estimated Vd using normal tissue assumption, error percent and standard deviation
% 51-54. Larsson filter:  original MTT,        estimated Vd, error percent and standard deviation
% 55-58. Larsson filter:  original MTT,        estimated MTT using normal tissue assumption, error percent and standard deviation
% 59-62. Larsson filter:  original E,          estimated E, error percent and standard deviation
% 63-66. Gaussian filter: original AIF delay,  estimated AIF delay using Gaussian, error percent and standard deviation
% 67-70. Sourbron method: original Flow,       estimated Flow, error percent and standard deviation 
% 71-74. Sourbron method: original Ki,         estimated Ki using 2CXM, error percent and standard deviation
% 75-78. Sourbron method: original Vb,         estimated Vb using 2CXM, error percent and standard deviation
% 79-82. Sourbron method: original Ve,         estimated Ve using 2CXM, error percent and standard deviation
Sim_Struct.num_results_parameters = 82;
Sim_Struct.results                = zeros(Sim_Struct.num_results_parameters,Sim_Struct.num_iterations);

%% Check parameters conflicts

% Ignoring gaussian calculation but wanting to iterate over it's parameters
if (Sim_Struct.Ignore_Gaussian_Calculation && (Sim_Struct.iterate_gaussian_sigma || Sim_Struct.iterate_gaussian_time_delay || Sim_Struct.iterate_gaussian_time_delay) )
    error('Set to ignore Gaussian calculation. Cant iterate over its parameters!');
end

% One-hot option only for de-convolution method
if (Sim_Struct.Use_Cyclic_Conv_4_ht_est + Sim_Struct.Use_Upsampling_Delay_Comp + Sim_Struct.Correct_estimation_due_to_delay > 1)
    error('More than 1 option for deconvolution is set!');
end

% Usampling + Cyclic can only be set in case Cyclic De-convolve is set
if (Sim_Struct.Use_Upsampling_and_Cyclic && ~Sim_Struct.Use_Cyclic_Conv_4_ht_est)
    error('Cyclic De-convolve is unset! Cant up-sample!');
end

% Use AIF delay corrcetion only in case we allowed it
if (Sim_Struct.Simple_AIF_Delay_Correct && ~Sim_Struct.Correct_estimation_due_to_delay)
    error('Cant correct for AIF delay, because correction is not enabled!');
end

if strcmp(Verbosity,'Full')
    display('-I- Finished Setting simulation parameters...');
end

end