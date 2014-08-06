classdef SimParams
    
    %Create params class with all needed algorithm parameters
    properties (Constant)
        
        %% AIF parameters - Parker's AIF
        
        % Factor for AIF to have realistic values like we see in MRI images
        r_factor = 0.0116;
        
        % AIF average population parameters
        A1    = 0.809;
        A2    = 0.330;
        
        % Changed to shift the boluses to the right
        T1    = 0.17046;
        T2    = 0.365;
        %T1    = 0.27046;
        %T2    = 0.465;
        
        sig1  = 0.0563;
        sig2  = 0.132;
        alpha = 1.050;
        beta  = 0.1685;
        s     = 38.078;
        tau   = 0.483;
        % Additional AIF parameters
        additional_AIF_delay_sec = 0; % Delay added to AIF before filtering
        additional_AIF_delay_min = SimParams.additional_AIF_delay_sec/60; % Translate to minutes
        
        %%
        
        % Gaussian filter parameters
        sigma          = 2/60;   % In minutes
        t_d            = 30/60;  % Time shift in minutes
        amplitude      = 1;      % Amplitude of Gaussian
        
        % Larsson filter parameters
        F_base     = 10; % When iterating
        F_single   = 60; % When a single iteration
        Vb_base    = 6;
        Vb_single  = 6;
        E_base     = 0.5;
        E_single   = 0.5;
        Ve_base    = 0.1;
        Ve_single  = 0.1;
        Hct_single = 0.38;
        
        % Number of FFT points
        num_FFT_points = 2000;
        
        % Time interval between samples
        sec_interval   = 2;                %[sec]
        Fs             = 1 / SimParams.sec_interval; %[Hz] - Sampling rate
        min_interval   = SimParams.sec_interval/60;  %[min]
        
        % Total simulation time
        total_sim_time_min  = 4; %[min]
        num_time_stamps     = round(SimParams.total_sim_time_min / SimParams.min_interval);
        
        % Time vector for AIF and Ct(t)
        time_vec_minutes    = (0:SimParams.num_time_stamps-1).* SimParams.min_interval;
        
        % Non-linear square fitting parameters
        FMS_TolFun          = 1e-11;
        FMS_MaxFunEvals     = 10000;
        FMS_MaxIter         = 10000;
        %FMS_Algorithm      = 'levenberg-marquardt'; % Does not supper lower-upper bounds
        FMS_Algorithm       = 'trust-region-reflective';
        algorithm_options   = optimset('TolFun',SimParams.FMS_TolFun,'MaxFunEvals',SimParams.FMS_MaxFunEvals,'MaxIter',SimParams.FMS_MaxIter,...
            'Display','off','Algorithm',SimParams.FMS_Algorithm);
        
        % Initial Guess for non-linear curve fitting for Gaussian (td=1sec, var=1[sec], amplitude = 1)
        Init_Guess_Gaussian = [1/60 1/60 0.5];
        
        % Set lower and upper bounds for parameters for non-linear fitting
        % Gaussian
        
        % Time delay -> 0 to 60 seconds
        % Sigma      -> 0 to 30 seconds
        % Amplitude  -> 0 to 3
        LowerBound_Gauss  = [0     0     0];
        UpperBound_Gauss  = [60/60 30/60 3];
        
        % Larsson parameters boundaries for non-linear curve fitting
        % Vb -> 0 to 100 [ml/100g]
        % E  -> 0 to 1
        % Ve -> 0 to 100 [ml/100g]
        LowerBound_Larsson  = [0   0     0];
        UpperBound_Larsson  = [100 1     100];
        
        % Determines SNR ( noise_var = mean(signal)/SNR_base )
        SNR_base = 15;
        
        % Additional parameters for simulation iterations
        Ignore_Gaussian_Calculation   = 1;      % Ignores all calculations relating to gaussian function
        Iterate_Murase_Tofts          = 0;      % Check Murase calculation for Tofts parameters
        Iterate_Murase_Tofts_num_iter = 100000; % number of iterations for Murase Tofts check
        One_Iteration_Murase_Tofts    = 0;
        
        % Drive differential equation
        Drive_Diff_Eq                 = 0;
        
        % Set if derivative matrix will be divided by time
        Derivative_Time_Devision = false;
        
        % Choose which Patlak estimation to take
        % Possible - 1. "Specified Points" 2. "All Points" 3. "Weighted Points"
        Patlak_Est = 'Specified Points';
        
        % Choose which filter estimation to use when calculating parameters (non-linear method)
        % Possible: 'Wiener', 'Ridge', 'Spline', 'Spline_1st', 'Spline_2nd'
        Filter_Est_Chosen = 'Spline_2nd';
        
        % Polynomial degree for basis splines
        poly_deg = 4;
        
        % Choose knots for splines (currently takes every 1 out of 2 points)
        knot_interval = 2;
        knots         = SimParams.time_vec_minutes(1:SimParams.knot_interval:end);
        
        % Lambda for regularization
        % 1. Ridge regression  2. Spline-no deriv
        % 3. Spline-1st deriv  4. Spline-2nd deriv
        %lambda_vec              = [9 132 0.0391 0.0141];
        lambda_vec_gauss        = [34.7 11 0.1 16.8];
        %lambda_vec              = [9 132 0.0391 0.4];
        lambda_vec_larss        = [4.7 0.9 0.1 0.4];
        %lambda_vec              = [9 132 0.0391 1.9];
        %lambda_vec              = [9 132 0.0391 0.0];
        %lambda_vec              = [9 132 0.0391 0.012]; 74.8
        normalize               = 1;  % Normalize flag for ridge() function
        
        % Choose whether to plot L curve
        plot_L_Curve = 0;
        
        
        
    end
    
    properties (Dependent)
        
        % Initiate plot parameters
        plot_flag = 0;
        plot_error_results_flag = 0;
        
        % Set number of iterations for simulation
        num_iterations = 15;
        %num_iterations = 1;
        
        % Do each iteration a few time and average results for better statistic information
        num_averages  = 2;
        %num_averages  = 1;
        
        
        % Determine according to what parameter to check simulations
        % Allow different than zero only in case we use more than 1 iteration
        iterate_SNR                   = 0 && (SimParams.num_iterations > 1);
        iterate_sec_interval          = 0 && (SimParams.num_iterations > 1);
        iterate_gaussian_sigma        = 0 && (SimParams.num_iterations > 1);
        iterate_gaussian_time_delay   = 0 && (SimParams.num_iterations > 1);
        iterate_gaussian_amplitude    = 0 && (SimParams.num_iterations > 1);
        iterate_F_larsson             = 1 && (SimParams.num_iterations > 1);
        iterate_Vb_larsson            = 0 && (SimParams.num_iterations > 1);
        iterate_E_larsson             = 0 && (SimParams.num_iterations > 1);
        iterate_Ve_larsson            = 0 && (SimParams.num_iterations > 1);
        
        % Initiate results matrix which will hold:
        % 1.     SNR ratio
        % 2.     Seconds interval of simulation
        % 3-6.   Gaussian:       original sigma, estimated sigma, error percent and standard deviation
        % 7-10.  Gaussian:       original delay time, estimated delay time, error percent and standard deviation
        % 11-14. Gaussian:       original amplitude, estimated amplitude, error percent and standard deviation
        % 15-18. Larsson filter: original Flow, estimated Flow, error percent and standard deviation
        % 19-22. AIF:            original delay, estimated delay, error percent and standard deviation
        % 23-26. Larsson filter: original Ki,  estimated Ki using Patlak, error percent and standard deviation
        % 27-30. Larsson filter: original Ki,  estimated Ki using 2CXM, error percent and standard deviation
        % 31-34. Larsson filter: original PS,  estimated PS, error percent and standard deviation
        % 35-38. Larsson filter: original Vb,  estimated Vb using Patlak, error percent and standard deviation
        % 39-42. Larsson filter: original Vb,  estimated Vb using 2CXM, error percent and standard deviation
        % 43-46. Larsson filter: original Vd,  estimated Vd, error percent and standard deviation
        % 47-50. Larsson filter: original Vd,  estimated Vd using normal tissue assumption, error percent and standard deviation
        % 51-54. Larsson filter: original MTT, estimated Vd, error percent and standard deviation
        % 55-58. Larsson filter: original MTT, estimated MTT using normal tissue assumption, error percent and standard deviation
        results       = zeros(58,SimParams.num_iterations);
        
    end
    
    
    methods
        
        function obj=update_params(obj)
            
            % Plot graphs when not iterating
            % When iterating, plot the error as a function of parameters
            if (obj.num_iterations == 1 && obj.num_averages == 1)
                obj.plot_flag = 1;
                obj.plot_error_results_flag = 0;
            else
                obj.plot_flag = 0;
                obj.plot_error_results_flag = 1;
            end
            
            
            % If checking Murase vs. Tofts, avoid averages iterations
            if (obj.Iterate_Murase_Tofts)
                obj.num_averages = 1;
            end
            
            % In case we want one iteration for graph display, change # to 1
            if (obj.One_Iteration_Murase_Tofts)
                obj.num_iterations = 1;
                obj.num_averages   = 1;
            end
            
            %% Check parameters conflicts
            % Ignoring gaussian calculation but wanting to iterate over it's parameters
            
            if (obj.Ignore_Gaussian_Calculation && (obj.iterate_gaussian_sigma || ...
                    obj.iterate_gaussian_time_delay || ...
                    obj.iterate_gaussian_time_delay) )
                error('-E- Set to ignore Gaussian calculation. Cant iterate over its parameters!');
            end
        end
        
    end
    
end