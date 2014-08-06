function [ Flow_vec, Delay_sec_vec, Delay_BiExp_Fit, t_delay_seconds_vec, sigma_seconds_vec, Amp_vec, Est_ht, ...
    calculated_gaussian, conv_result_ht, conv_result_gaussian, RMS_ht, RMS_gauss, RMS_params,...
    calculated_double_gaussian, conv_result_double_gaussian, double_gaussian_param_vec, ...
    RMS_double_gauss, RMS_params_double_gauss, Ki_vec, Vb_vec, Ve_vec, MTT_vec, Ki_Patlak_vec, Vb_Patlak_vec, MTT_Patlak_vec ] ...
    = Get_Ht_Deconvolving( AIF, Ct , sec_interval, Output_directory, Subject_name, FORCE)

%Get_Ht_Deconvolving Extracting possible h_t filter
%   The function gets AIF(nT), Ct(num_voxels,nT) - the arterial input function and the
%   tissue conventration (in minutes metric).
%   In addition it gets the time interval between samples (sec_interval).
%   It uses deconvolution (Wiener Filter / Tichonov Regularization) to estimate the convoultion
%   filter ( h(t) ).
%   The SNR approximation is based on experimental results.
%   After h(t) aprox. , it fits a gaussian/larss filter using non-linear fitting method.
%   The output consist of:
%   - t_delay, sigma and Amp, the gaussian parameters
%   - h_t, the deconvolved filter estimation
%   - calculated_gaussian, the fitted gaussian to ht
%   - conv_result_ht, the result of AIF convolved with h_t
%   - conv_result_gaussian, the result of AIF convolved with gaussian
%   - RMS_ht, the root mean square error of conv_result with h_t and Ct.
%   - RMS_gauss, the root mean square error of conv_result with estimated gaussian and Ct.
%   - RMS_params, the root mean square error of the gaussian using the parameters
%          and h(t)  (if too big, h(t) is not really a gaussian)
%   - calculated_double_gaussian, the fitted double gaussian to ht
%   - conv_result_double_gaussian,  the result of AIF convolved with double gaussian
%   - double_gaussian_param_vec, all the double gaussian parameters (means, vars and amplitudes).
%   - RMS_double_gauss, the root mean square error of conv_result with estimated double gaussian and Ct.
%   - RMS_params_2, the root mean square error of the double gaussian using the parameters
%                   and h(t)  (if too big, h(t) is not really double gaussian)
%   - Ki   - estimation for permeability according to Larsson
%   - Vb   - estimation for blood volume according to Larsson
%   - Ve   - estimation for EES volume according to Larsson
%   - MTT  - estimation for Mean Transient Time according to Larsson
%   - Ki_Patlak_vec  - estimation for permeability according to Patlak
%   - Vb_Patlak_vec  - estimation for blood volume according to Patlak
%   - MTT_Patlak_vec - estimation for Mean Transient Time according to Patlak


%% Initiate parameters

% Choose what to use/not use
USE_ONE_GAUSSIAN    = false;
USE_DOUBLE_GAUSSIAN = false;
USE_WIENER          = false;
USE_TICHONOV        = true;
% Try to correct for possible AIF delay
Delay_Correct       = false;

% Time interval between samples
Fs             = 1 / sec_interval; %[Hz] - Sampling rate
min_interval   = sec_interval/60;  %[min]

% Get number of time stamps (vector length)
num_time_stamps    = size(Ct,2);
num_voxels         = size(Ct,1);

% Time vector for AIF and Ct(t)
time_vec_minutes   = double((1:num_time_stamps).* min_interval);
time_vec_minutes_T = transpose(time_vec_minutes);

% Non-linear square fitting parameters
FMS_TolFun        = 1e-11;
FMS_MaxFunEvals   = 10000;
FMS_MaxIter       = 10000;
algorithm_options = optimset('TolFun',FMS_TolFun,'MaxFunEvals',FMS_MaxFunEvals,'MaxIter',FMS_MaxIter,'Display','off');

% Initial Guess (td=1sec, var=1, amp=0)
init_guess   = [1/60 1/60 0];
init_guess_2 = [1/60 1/60 0 3/60 1/60 0];

% Set lower and upper bounds for parameters
% Time delay -> 0 to 3 seconds
% Sigma      -> 0 to 1 seconds
% Amplitude  -> 0 to 1
LowerBound  = [0     0    0];
UpperBound  = [30/60 1/60 1];

% For double gaussian
LowerBound_2  = [0     0    0  0    0     0];
UpperBound_2  = [30/60 1/60 1 4*60/60 60/60 1];


% Larsson parameters boundaries for non-linear curve fitting
% Vb -> 0 to 100 [ml/100g]
% E  -> 0 to 1
% Ve -> 0 to 100 [ml/100g]
LowerBound_Larsson  = [0   0     0];
UpperBound_Larsson  = [100 1     100];

%% Estimating h(t) by Wiener filter / Tikhonov Regularization
Est_ht          = zeros(num_voxels,num_time_stamps);
Delay_BiExp_Fit = zeros(1,num_voxels);

% -------------- Wiener ---------------------------

if (USE_WIENER)
    [Est_ht] = Wiener_Filter( min_interval*AIF, Ct, Fs);
end

% -------------- Tikhonov ---------------------------

if (USE_TICHONOV)
    
    % Parameters
    %lambda_vec              = [9 132 0.0391 0.0141];
    lambda_vec              = [4.7 0.9 0.1 0.2];
    normalize               = 1;  % Normalize flag for ridge() function
    % Choose knots for splines (currently takes every 1 out of 2 points)
    knot_interval           = 2;
    knots                   = time_vec_minutes(1:knot_interval:end);
    % Create the the spline basis matrix
    poly_deg                = 4;
    display('-I- Create B-splines matrix...');
    B_mat                   = Create_B_matrix(knots,time_vec_minutes,poly_deg-1);
    
    % Create convolution indices
    [ Conv_X ] = Convolution_Matrix( min_interval, AIF );
    
    % Check if already calculated Ht
    Mat_File_Ht = [Output_directory 'Estimated_Tichonov_Ht_' Subject_name '.mat'];
    
    if(exist(Mat_File_Ht,'file') && ~FORCE)
        load(Mat_File_Ht);
    else
        
        display('--------------------------------------------------------');
        display('-I- Starting h(t) estimation using regularization...');
        display('--------------------------------------------------------');
        
        %for j=1:num_voxels
        parfor j=1:num_voxels
            
            % Deconvolution by regularization for larsson's filter
            [Est_ht(j,:), Delay_BiExp_Fit(j)] = Regularization_Methods(Ct(j,:)', AIF', Conv_X, time_vec_minutes, lambda_vec, normalize, min_interval, B_mat, Delay_Correct );
            
            % Report after each 1000 voxels
            if ( mod(j,1000) == 0 )
                
                display(sprintf('Finished Regularization_Methods for %d voxels...',j));
                
                remaining_voxels = num_voxels - j;
                
                fprintf('Number of remaining voxels: %d .\n',remaining_voxels);
                
            end
            
        end
        save(Mat_File_Ht,'Est_ht');
        
    end
    
end

% Create the transpose for specific functions
Est_ht_T = transpose(Est_ht);

% The analytic funcational of a gaussian function (x(3) is the amplitude)
Gaussian_function        = @(x,t) Gaussian      ( t,x(1),x(2),x(3)                );
Double_Gaussian_function = @(x,t) DoubleGaussian( t,x(1),x(2),x(3),x(4),x(5),x(6) );
% Hematocrit according to Larsson's article
Hct                      = 0.38;

% Gaussian parameters
t_delay_in_min_vec            = zeros(1,num_voxels);
t_delay_seconds_vec           = zeros(1,num_voxels);
est_var_vec                   = zeros(1,num_voxels);
sigma_in_min_vec              = zeros(1,num_voxels);
sigma_seconds_vec             = zeros(1,num_voxels);
Amp_vec                       = zeros(1,num_voxels);
double_gaussian_param_vec     = zeros(6,num_voxels);
calculated_gaussian           = zeros(num_voxels,num_time_stamps);
calculated_double_gaussian    = zeros(num_voxels,num_time_stamps);

% Larsson parameters
Flow_vec                      = zeros(1,num_voxels);
Delay_sec_vec                 = zeros(1,num_voxels);
Ki_vec                        = zeros(1,num_voxels);
Vb_vec                        = zeros(1,num_voxels);
Ve_vec                        = zeros(1,num_voxels);
MTT_vec                       = zeros(1,num_voxels);
Ki_Patlak_vec                 = zeros(1,num_voxels);
Vb_Patlak_vec                 = zeros(1,num_voxels);
MTT_Patlak_vec                = zeros(1,num_voxels);


%% Estimate parameters for all voxels by curve fitting

% start_time = now;
% tic;
% sub_iter_num = 0;
% duration_of_100_average = 0;

% Define par struct for performance analysis
%par_iter = Par(num_voxels);

% Check if already estimated parameters for all voxels
Mat_File_Perfusion_Parameters = [Output_directory 'Estimated_Perfusion_Parameters_' Subject_name '.mat'];

if(exist(Mat_File_Perfusion_Parameters,'file') && ~FORCE)
    load(Mat_File_Perfusion_Parameters);
else
    
    display('--------------------------------------------------------');
    display('-I- Starting non linear parameters estimation...');
    display('--------------------------------------------------------');
    
    % Add needed functions to parallel workers
    myPool = gcp;
    myFiles = {'Gaussian.m', 'DoubleGaussian.m'};
    addAttachedFiles(myPool, myFiles);
    
    parfor j=1:num_voxels
    %for j=1:num_voxels
        
        %         display(sprintf('-I- Voxel number: %d',j));
        %         if (j==48)
        %             a = 1;
        %         end
        
        
        % Start timing Calculation
        %Par.tic;
        

        
        % lsqcurvefit parameters are:
        % analytic function, initial parameters, time vector, data points ,lower
        % and upper bounds and algorithm options
        
        if (USE_ONE_GAUSSIAN)
            [est_params,residue_norm,residual,exitflag,algo_info] = ...
                lsqcurvefit(Gaussian_function,init_guess,time_vec_minutes_T,Est_ht_T(:,j),LowerBound,UpperBound,algorithm_options);
        end
        
        
        if (USE_DOUBLE_GAUSSIAN)
            [est_params_2,residue_norm_2,residual_2,exitflag_2,algo_info_2] = ...
                lsqcurvefit(Double_Gaussian_function,init_guess_2,time_vec_minutes_T,Est_ht_T(:,j),LowerBound_2,UpperBound_2,algorithm_options);
        end
        
        est_params                    = zeros(1,3);
        est_params_2                  = zeros(1,6);
        
        % Put parameters in wanted variables
        t_delay_in_min_vec(j)  = est_params(1);
        t_delay_seconds_vec(j) = 60 * t_delay_in_min_vec(j); % Convert delay time to seconds
        est_var_vec(j)         = est_params(2);
        sigma_in_min_vec(j)    = sqrt(est_var_vec(j));
        sigma_seconds_vec(j)   = 60 * sigma_in_min_vec(j); % Convert sigma to seconds
        Amp_vec(j)             = est_params(3);
        
        % Larsson parameters
        Flow_vec(j)      = max(Est_ht_T(:,j));
        
        %% ----------------------- PATLAK --------------------------------------
        
        %Use patlak to get initial parameters estimation
        
        
        % Ignore points where Ca(t) is very small (because then the result is unstable)
        Y_vec_Vb         = Ct(j,:) ./ AIF;                        %[mL/100g]
        X_vec            = cumtrapz(time_vec_minutes,AIF) ./ AIF; %[min]
        [val, bolus_idx] = max(diff(AIF));
        base_value       = mean(AIF(1:bolus_idx-1));
        mult_val_Thresh  = 3;
        Threshold        = mult_val_Thresh*base_value;
        stable_idx       = find(AIF > Threshold);
        
        % Take the stable points out of the vector
        X_vec            = X_vec(stable_idx);
        Y_vec_Vb         = Y_vec_Vb(stable_idx);
        
        % Remove Zeros/NaNs/Infs caused because of division by 0
        nan_indices           = find(isnan(Y_vec_Vb));
        inf_indices           = find(~isfinite(Y_vec_Vb));
        Y_vec_Vb(nan_indices) = [];
        Y_vec_Vb(inf_indices) = [];
        X_vec(nan_indices)    = [];
        X_vec(inf_indices)    = [];
        
        % Fine straight line coefficent
        [linear_params]       = polyfit(X_vec,Y_vec_Vb,1);
        
        est_Ki_Patlak_noise   = linear_params(1); %a in ax+b
        est_Vb_Patlak_noise   = linear_params(2); %b in ax+b
        
        % Handle zero value
        if ( Flow_vec(j) > 0 )
            E_Patlak_est          = est_Ki_Patlak_noise / Flow_vec(j); % E = Ki / F
            est_MTT_Patlak        = est_Vb_Patlak_noise / Flow_vec(j);
        else
            est_Vb_Patlak_noise   = 0;
            E_Patlak_est          = 0;
            est_MTT_Patlak        = 0;
        end
        
        
        %% ----------------------- PATLAK END --------------------------------------
        
        % Initial Guess for non-linear curve fitting for Larsson (Vb, E, Ve)
        Init_Guess_Larsson = double( [est_Vb_Patlak_noise E_Patlak_est 10] );
        
        est_F = double(Flow_vec(j));
        % The analytic funcational of a Larsson function
        Larsson_function         = @(x,t) Larsson_Filter( t, est_F, x(1), x(2), x(3), Hct);
        
        % Calculate parameters only in case F is not 0 (or else there is a
        % problem)
        
        if (est_F~=0)
            [est_params_Larsson_noise,residue_norm_Larsson_noise,residual_Larsson_noise,exitflag_Larsson_noise,algo_info_Larsson_noise] = ...
                lsqcurvefit(Larsson_function,Init_Guess_Larsson,time_vec_minutes_T,Est_ht_T(:,j)/Flow_vec(j),...
                LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
            
            % Assigning two compartment parameters estimation
            Ki_vec(j)         = est_params_Larsson_noise(2)* Flow_vec(j);
            Vb_vec(j)         = est_params_Larsson_noise(1);
            Ve_vec(j)         = est_params_Larsson_noise(3);
            % Estimate MTT
            est_IRF           = Est_ht_T(:,j) / Flow_vec(j);
            est_MTT_noise     = cumtrapz(time_vec_minutes,est_IRF);
            MTT_vec(j)        = est_MTT_noise(end);
            
            % Assigning Patlak parameters estimation
            Ki_Patlak_vec(j)  = est_Ki_Patlak_noise;
            Vb_Patlak_vec(j)  = est_Vb_Patlak_noise;
            MTT_Patlak_vec(j) = est_MTT_Patlak;
            
            % Delay of the AIF will be calculated according to the place of the maximum value of F*IRF
            max_index        = find( Est_ht_T(:,j) == Flow_vec(j) );
            % Translate to minutes
            Delay_sec_vec(j) =  (max_index - 1) * min_interval * 60;
            
            % Convert sigma to variance and minutes to seconds
            delay_1_seconds = est_params_2(1) * 60;
            sigma_1_seconds = sqrt(est_params_2(2)) * 60;
            amp_1           = est_params_2(3);
            delay_2_seconds = est_params_2(4) * 60;
            sigma_2_seconds = sqrt(est_params_2(5)) * 60;
            amp_2           =  est_params_2(6);
            
            double_gaussian_param_vec(:,j) = [delay_1_seconds sigma_1_seconds amp_1 delay_2_seconds sigma_2_seconds amp_2];
            
            
            % Calculate gaussian out of parameters
            calculated_gaussian(j,:)        = Gaussian(time_vec_minutes, t_delay_in_min_vec(j), sigma_in_min_vec(j)^2, Amp_vec(j));
            
            % Calculate double gaussian out of parameters
            calculated_double_gaussian(j,:) = DoubleGaussian(time_vec_minutes, est_params_2(1), est_params_2(2), est_params_2(3) ...
                ,  est_params_2(4), est_params_2(5), est_params_2(6));
            
        end
        
        % Report after each 100 voxels
        if ( mod(j,100) == 0 )
            
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
    
    
    save(Mat_File_Perfusion_Parameters,'Flow_vec','t_delay_in_min_vec','t_delay_seconds_vec','est_var_vec','sigma_in_min_vec','sigma_seconds_vec',...
        'Amp_vec','Ki_vec','Vb_vec','Ve_vec','MTT_vec','Ki_Patlak_vec','Vb_Patlak_vec',...
        'MTT_Patlak_vec','Delay_sec_vec','double_gaussian_param_vec','calculated_gaussian','calculated_double_gaussian');
    
end


% Display performance analysis
% stop(par_iter);
% figure;
% plot(par_iter);

% Calculate RMS of convolution result comparing to Ct(t)
RMS_params              = sqrt( sum( (calculated_gaussian        - Est_ht).^2 ,2) );
RMS_params_double_gauss = sqrt( sum( (calculated_double_gaussian - Est_ht).^2 ,2) );

%% Filter AIF through kernel

% Filter the AIF with the estimated ht
%conv_result_ht       = filter(Est_ht,1,AIF);
conv_result_ht       = filter(AIF*min_interval,1,Est_ht,[],2);

% Filter the AIF with the gaussian kernel
%conv_result_gaussian = filter(calculated_gaussian*min_interval,1,AIF);
conv_result_gaussian        = filter(AIF*min_interval,1,calculated_gaussian,[],2);
conv_result_double_gaussian = filter(AIF*min_interval,1,calculated_double_gaussian,[],2);


% Zero negative values
conv_result_ht(conv_result_ht<0)             = 0;
conv_result_gaussian(conv_result_gaussian<0) = 0;
conv_result_double_gaussian(conv_result_double_gaussian<0) = 0;

% Calculate RMS of convolution results comparing to Ct(t)
RMS_ht           = sqrt( sum( (Ct - conv_result_ht).^2 , 2) );
RMS_gauss        = sqrt( sum( (Ct - conv_result_gaussian).^2, 2) );
RMS_double_gauss = sqrt( sum( (Ct - conv_result_double_gaussian).^2, 2) );

end

