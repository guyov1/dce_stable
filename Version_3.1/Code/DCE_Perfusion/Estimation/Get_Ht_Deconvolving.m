function [ returnStruct ] = Get_Ht_Deconvolving(Sim_Struct, AIF, Ct , Output_directory, Subject_name, Force_RealData_Calc, Verbosity)

%Get_Ht_Deconvolving Extracting possible h_t filter
%   The function gets AIF(nT), Ct(num_voxels,nT) - the arterial input function and the tissue conventration (in minutes metric).
%   In addition it gets the time interval between samples (sec_interval).
%   It uses deconvolution (Wiener Filter / Tichonov Regularization) to estimate the convoultion filter ( h(t) ).
%   The SNR approximation is based on experimental results.
%   After h(t) aprox. , it fits a gaussian/larss filter using non-linear fitting method.
%
%   The output consist of:
%   - t_delay, sigma and Amp, the gaussian parameters
%   - h_t, the deconvolved filter estimation
%   - calculated_gaussian, the fitted gaussian to ht
%   - conv_result_ht, the result of AIF convolved with h_t
%   - conv_result_gaussian, the result of AIF convolved with gaussian
%   - RMS_ht, the root mean square error of conv_result with h_t and Ct.
%   - RMS_gauss, the root mean square error of conv_result with estimated gaussian and Ct.
%   - RMS_params_Gauss, the root mean square error of the gaussian using the parameters
%          and h(t)  (if too big, h(t) is not really a gaussian)
%   - calculated_double_gaussian, the fitted double gaussian to ht
%   - conv_result_double_gaussian,  the result of AIF convolved with double gaussian
%   - double_gauss_param, all the double gaussian parameters (means, vars and amplitudes).
%   - RMS_double_gauss, the root mean square error of conv_result with estimated double gaussian and Ct.
%   - RMS_params_double_gauss, the root mean square error of the double gaussian using the parameters
%                   and h(t)  (if too big, h(t) is not really double gaussian)
%   - Ktrans   - estimation for permeability according to Larsson
%   - Vb   - estimation for blood volume according to Larsson
%   - Ve   - estimation for EES volume according to Larsson
%   - MTT  - estimation for Mean Transient Time according to Larsson
%   - Ktrans_Patlak  - estimation for permeability according to Patlak
%   - Vb_Patlak      - estimation for blood volume according to Patlak
%   - MTT_Patlak     - estimation for Mean Transient Time according to Patlak


%% Initiate parameters

Adjusted_Larsson_Model          = Sim_Struct.Adjusted_Larsson_Model;
Patlak_Est_Type                 = Sim_Struct.Patlak_Est_Type;
Parallel_Real_Data_Est          = Sim_Struct.Parallel_Real_Data_Est;
time_vec_minutes                = Sim_Struct.time_vec_minutes;
min_interval                    = Sim_Struct.min_interval;
num_time_stamps                 = Sim_Struct.num_time_stamps;
num_voxels                      = Sim_Struct.num_voxels;
lambda_vec_larss                = Sim_Struct.lambda_vec_larss;
normalize                       = Sim_Struct.normalize;  % Normalize flag for ridge() function
plot_L_Curve                    = Sim_Struct.plot_L_Curve;
Derivative_Time_Devision        = Sim_Struct.Derivative_Time_Devision;
Upsampling_resolution           = Sim_Struct.Upsampling_resolution;
Max_Time_Delay                  = Sim_Struct.Max_Time_Delay;
Min_Time_Delay                  = Sim_Struct.Min_Time_Delay;
RMS_Smooth_Around_Bolus         = Sim_Struct.RMS_Smooth_Around_Bolus;
RMS_Smooth                      = Sim_Struct.RMS_Smooth;
Correct_estimation_due_to_delay = Sim_Struct.Correct_estimation_due_to_delay;
Simple_AIF_Delay_Correct        = Sim_Struct.Simple_AIF_Delay_Correct;
LowerBound_Larsson              = Sim_Struct.LowerBound_Larsson;
UpperBound_Larsson              = Sim_Struct.UpperBound_Larsson;
Hct                             = Sim_Struct.Hct_single;
Diff_From_Bolus                 = Sim_Struct.Diff_From_Bolus;        % The difference in seconds from the bolus to look on
Vb_low                          = Sim_Struct.Vb_low;
algorithm_options               = Sim_Struct.algorithm_options;
Use_Upsampling_Delay_Comp       = Sim_Struct.Use_Upsampling_Delay_Comp;
BiExp2CTC_RMS_Ratio             = Sim_Struct.BiExp2CTC_RMS_Ratio;
Filter_Est_Chosen               = Sim_Struct.Filter_Est_Chosen;
USE_WIENER                      = Sim_Struct.USE_WIENER;
USE_TICHONOV                    = Sim_Struct.USE_TICHONOV;
RealData_Flag                   = Sim_Struct.RealData_Flag;
plot_flag                       = Sim_Struct.plot_flag;
filter_type                     = Sim_Struct.filter_type;

% Put dummy values that are used in simulation
Sim_Ct_T                 = NaN;
Conv_Matrix_no_noise     = NaN;
B_PCA                    = NaN;
idx_fig                  = 1;
% filter_type              = 'Larss';


%% Estimating h(t) by Wiener filter / Tikhonov Regularization
Est_IRF_with_Delay              = zeros(num_voxels,num_time_stamps);
Est_IRF_no_Delay                = zeros(num_voxels,num_time_stamps);
est_delay_by_AIF_correct        = zeros(1,num_voxels);

if USE_WIENER
    [Est_IRF_no_Delay] = Wiener_Filter( min_interval*AIF, Ct, Fs);
elseif USE_TICHONOV
    
    % Choose knots for splines (currently takes every 1 out of 2 points)
    knot_interval           = Sim_Struct.knot_interval;
    knots                   = time_vec_minutes(1:knot_interval:end);
    % Create the the spline basis matrix
    poly_deg                = Sim_Struct.poly_deg;
    B_mat                   = Create_B_matrix(knots,time_vec_minutes,poly_deg-1);
    % Create convolution indices
    [ Conv_Matrix ] = Convolution_Matrix( min_interval, AIF );
    
    % Check if already calculated Ht
    Mat_File_Ht = [Output_directory 'Estimated_Tichonov_Ht_' Subject_name '.mat'];
    
    % --------------- AIF delay correction parameters ------------------
    Delay_Correct_Struct                                    = struct;
    Delay_Correct_Struct.normalize                          = normalize;
    Delay_Correct_Struct.B_mat                              = B_mat;
    Delay_Correct_Struct.PCA_B_mat                          = NaN;
    Delay_Correct_Struct.plot_L_Curve                       = plot_L_Curve;
    Delay_Correct_Struct.Derivative_Time_Devision           = Derivative_Time_Devision;
    Delay_Correct_Struct.lambda_vec_larss                   = lambda_vec_larss;
    Delay_Correct_Struct.min_interval                       = min_interval;
    Delay_Correct_Struct.time_vec_minutes                   = time_vec_minutes;
    Delay_Correct_Struct.Upsampling_resolution              = Upsampling_resolution;
    Delay_Correct_Struct.Max_Time_Delay                     = Max_Time_Delay;
    Delay_Correct_Struct.Min_Time_Delay                     = Min_Time_Delay;
    Delay_Correct_Struct.Use_Upsampling_Delay_Comp          = Use_Upsampling_Delay_Comp;
    Delay_Correct_Struct.LowerBound_Larsson                 = LowerBound_Larsson;
    Delay_Correct_Struct.UpperBound_Larsson                 = UpperBound_Larsson;
    Delay_Correct_Struct.algorithm_options                  = algorithm_options;
    Delay_Correct_Struct.Hct                                = Hct; % Try to read it from patient
    Delay_Correct_Struct.RMS_Smooth_Around_Bolus            = RMS_Smooth_Around_Bolus;
    Delay_Correct_Struct.RMS_Smooth                         = RMS_Smooth;
    Delay_Correct_Struct.Diff_From_Bolus                    = Diff_From_Bolus;
    Delay_Correct_Struct.additional_AIF_delay_sec           = 0;
    Delay_Correct_Struct.BiExp2CTC_RMS_Ratio                = BiExp2CTC_RMS_Ratio;
    Delay_Correct_Struct.plot_flag                          = false;
    Delay_Correct_Struct.Adjusted_Larsson_Model             = Adjusted_Larsson_Model;
    Delay_Correct_Struct.Filter_Est_Chosen                  = Filter_Est_Chosen;
    Delay_Correct_Struct.Vb_low                             = Vb_low;
    Delay_Correct_Struct.RealData_Flag                      = RealData_Flag;
    Delay_Correct_Struct.Simple_AIF_Delay_Correct           = Simple_AIF_Delay_Correct;
    Delay_Correct_Struct.Patlak_Est_Type                    = Patlak_Est_Type;
    Delay_Correct_Struct.Ktrans                             = NaN;       % Simulation ground truth values
    Delay_Correct_Struct.Vb_larss                           = NaN; % Simulation ground truth values
    Delay_Correct_Struct.init_Ve_guess                      = Sim_Struct.init_Ve_guess;
    Delay_Correct_Struct.FMS_Algorithm                      = Sim_Struct.FMS_Algorithm;
    Delay_Correct_Struct.LQ_Model_AIF_Delay_Correct         = Sim_Struct.LQ_Model_AIF_Delay_Correct;
    Delay_Correct_Struct.filter_type                        = Sim_Struct.filter_type;
    
    Delay_Correct_Ht_Struct                                 = struct;
    Delay_Correct_Ht_Struct.Sim_AIF_with_noise_Regul        = AIF;
    Delay_Correct_Ht_Struct.Sim_Ct_larss_Regul              = NaN;
    Delay_Correct_Ht_Struct.Conv_X_no_noise                 = NaN;
    
    % The iteration and average number of specific constant parameters
    iter_num                                                = 1;
    avg_num                                                 = 1;
    % ------------------------------------------------------------
    
    if(exist(Mat_File_Ht,'file') && ~Force_RealData_Calc)
        load(Mat_File_Ht);
        display('--------------------------------------------------------');
        display('-I- Loading pre-calculated h(t) estimation using regularization...');
        display('--------------------------------------------------------');
    else
        
        
        display('--------------------------------------------------------');
        display('-I- Starting h(t) estimation using regularization...');
        display('--------------------------------------------------------');
        
        AIF_delay_corrected = zeros(size(Ct));
        
        if Parallel_Real_Data_Est
            
            parfor j=1:num_voxels
                tic;
                
                [ ~, ~, ~, b_spline_result_2nd_deriv_no_Delay, ~, ~, ~, ~ ] =  ...
                    Regularization_Methods_Simulation( Sim_Ct_T, Ct(j,:)', Conv_Matrix, Conv_Matrix_no_noise, time_vec_minutes, lambda_vec_larss, normalize, min_interval, B_mat, B_PCA, plot_L_Curve, idx_fig, filter_type, Derivative_Time_Devision, plot_flag, RealData_Flag );
                
                Est_IRF_no_Delay(j,:) = b_spline_result_2nd_deriv_no_Delay;
                
                % Estimate delay and correct h(t) estimation if it seems we have delay in AIF
                if Correct_estimation_due_to_delay
                    [est_delay_by_AIF_correct(j), AIF_delay_corrected(j ,:), Est_IRF_with_Delay(j,:), ~] = ...
                        AIF_Delay_Correct(Delay_Correct_Struct, Delay_Correct_Ht_Struct, Est_IRF_no_Delay(j,:), Ct(j,:)', Verbosity, iter_num, avg_num, idx_fig);
                else
                    Est_IRF_with_Delay(j ,:)  = Est_IRF_no_Delay(j,:);
                    AIF_delay_corrected(j ,:) = AIF;
                end
                
                % if Use_Cyclic_Conv_4_ht_est
                % Cyclic_Deconvolve( Sim_Struct, Verbosity, iter_num, avg_num, idx_fig )
                
                % Report after each 1000 voxels
                if ( mod(j,1000) == 0 )
                    display(sprintf('Finished Regularization_Methods_Simulation for 1000 voxels in %d minutes...', toc/60));
                    remaining_voxels = num_voxels - j;
                    fprintf('Number of remaining voxels: %d .\n',remaining_voxels);
                end
                
            end
            
        else
            for j=1:num_voxels
                tic;
                
                [ ~, ~, ~, b_spline_result_2nd_deriv_no_Delay, ~, ~, ~, ~ ] =  ...
                    Regularization_Methods_Simulation( Sim_Ct_T, Ct(j,:)', Conv_Matrix, Conv_Matrix_no_noise, time_vec_minutes, lambda_vec_larss, normalize, min_interval, B_mat, B_PCA, plot_L_Curve, idx_fig, filter_type, Derivative_Time_Devision, plot_flag, RealData_Flag );
                
                Est_IRF_no_Delay(j,:) = b_spline_result_2nd_deriv_no_Delay;
                
                % Estimate delay and correct h(t) estimation if it seems we have delay in AIF
                if Correct_estimation_due_to_delay
                    [est_delay_by_AIF_correct(j), AIF_delay_corrected(j ,:), Est_IRF_with_Delay(j,:), ~] = ...
                        AIF_Delay_Correct(Delay_Correct_Struct, Delay_Correct_Ht_Struct, Est_IRF_no_Delay(j,:), Ct(j,:)', Verbosity, iter_num, avg_num, idx_fig);
                else
                    Est_IRF_with_Delay(j ,:)  = Est_IRF_no_Delay(j,:);
                    AIF_delay_corrected(j ,:) = AIF;
                end
                
                % if Use_Cyclic_Conv_4_ht_est
                % Cyclic_Deconvolve( Sim_Struct, Verbosity, iter_num, avg_num, idx_fig )
                
                % Report after each 1000 voxels
                if ( mod(j,1000) == 0 )
                    display(sprintf('Finished Regularization_Methods_Simulation for 1000 voxels in %d minutes...', toc/60));
                    remaining_voxels = num_voxels - j;
                    fprintf('Number of remaining voxels: %d .\n',remaining_voxels);
                end
                
            end
            
        end
        
        save(Mat_File_Ht, 'Est_IRF_with_Delay','Est_IRF_no_Delay', 'est_delay_by_AIF_correct', 'AIF_delay_corrected');
        
    end
    
end

%% Estimate parameters for all voxels by curve fitting

% Check if already estimated parameters for all voxels
Mat_File_Perfusion_Parameters = [Output_directory 'Estimated_Perfusion_Parameters_' Subject_name '.mat'];

if(exist(Mat_File_Perfusion_Parameters,'file') && ~Force_RealData_Calc)
    load(Mat_File_Perfusion_Parameters);
else
    [   Flow_with_Delay, Flow_no_Delay, Delay_sec_by_Max_Val_with_Delay, Delay_sec_by_Max_Val_no_Delay, ...
        fitted_larsson_with_Delay, fitted_larsson_no_Delay, fitted_larsson_with_Delay_High_F, fitted_larsson_no_Delay_High_F, fitted_larsson_with_Delay_no_Ve, fitted_larsson_no_Delay_no_Ve,...
        fitted_larsson_with_Delay_no_E, fitted_larsson_no_Delay_no_E, fitted_larsson_with_Delay_no_E_High_F, fitted_larsson_no_Delay_no_E_High_F, ...
        fitted_gaussian, fitted_double_gaussian, gaussian_param, double_gauss_param, ...
        Ktrans_with_Delay, Ktrans_no_Delay, ...
        E_with_Delay, E_no_Delay, E_with_Delay_no_Ve, E_no_Delay_no_Ve, Ktrans_with_Delay_High_F, Ktrans_no_Delay_High_F,  ...
        Vb_with_Delay, Vb_no_Delay, Vb_with_Delay_High_F, Vb_no_Delay_High_F, Vb_with_Delay_no_Ve, Vb_no_Delay_no_Ve,...
        Vb_with_Delay_no_E, Vb_no_Delay_no_E, Vb_with_Delay_no_E_High_F, Vb_no_Delay_no_E_High_F, ...
        Ve_with_Delay, Ve_no_Delay, Ve_with_Delay_High_F, Ve_no_Delay_High_F, MTT_with_Delay, MTT_no_Delay, ...
        Ktrans_Patlak_with_Delay, Ktrans_Patlak_no_Delay, Vb_Patlak_with_Delay, Vb_Patlak_no_Delay, ...
        MTT_Patlak_with_Delay, MTT_Patlak_no_Delay ] ...
        = nonLinParamEst(Sim_Struct, Est_IRF_with_Delay', Est_IRF_no_Delay', Ct, AIF_delay_corrected', AIF', idx_fig, Parallel_Real_Data_Est );
    
    save(Mat_File_Perfusion_Parameters,'Flow_with_Delay','Flow_no_Delay',...
        'Ktrans_with_Delay','Ktrans_no_Delay', 'E_with_Delay', 'E_no_Delay', 'E_with_Delay_no_Ve', 'E_no_Delay_no_Ve', 'Ktrans_with_Delay_High_F', 'Ktrans_no_Delay_High_F', ...
        'Vb_with_Delay', 'Vb_no_Delay', 'Vb_with_Delay_High_F', 'Vb_no_Delay_High_F', 'Vb_with_Delay_no_Ve', 'Vb_no_Delay_no_Ve', 'Vb_with_Delay_no_E', 'Vb_no_Delay_no_E', 'Vb_with_Delay_no_E_High_F', 'Vb_no_Delay_no_E_High_F', ...
        'Ve_with_Delay', 'Ve_no_Delay', 'Ve_with_Delay_High_F', 'Ve_no_Delay_High_F', 'MTT_with_Delay', 'MTT_no_Delay',...
        'Ktrans_Patlak_with_Delay','Ktrans_Patlak_no_Delay', 'Vb_Patlak_with_Delay','Vb_Patlak_no_Delay',...
        'MTT_Patlak_with_Delay','MTT_Patlak_no_Delay','Delay_sec_by_Max_Val_with_Delay',...
        'Delay_sec_by_Max_Val_no_Delay','gaussian_param', 'double_gauss_param','fitted_gaussian','fitted_double_gaussian', ...
        'fitted_larsson_with_Delay','fitted_larsson_no_Delay', 'fitted_larsson_with_Delay_High_F', 'fitted_larsson_no_Delay_High_F', 'fitted_larsson_with_Delay_no_Ve', 'fitted_larsson_no_Delay_no_Ve',  'fitted_larsson_no_Delay_no_E', 'fitted_larsson_with_Delay_no_E', 'fitted_larsson_with_Delay_no_E_High_F', 'fitted_larsson_no_Delay_no_E_High_F');
    
end

% Calculate RMS of convolution result comparing to Ct(t)
RMS_params_Gauss        = sqrt( sum( (fitted_gaussian        - Est_IRF_no_Delay).^2 ,2) );
RMS_params_double_gauss = sqrt( sum( (fitted_double_gaussian - Est_IRF_no_Delay).^2 ,2) );

%% Filter AIF through kernel

% Filter the AIF with the estimated ht - conv_result_ht = filter(Est_ht*min_interval,1,AIF)
if Sim_Struct.ignore_time_delta
    min_interval = 1;
end

% Fit without bi-exp fit
conv_result_no_Delay_IRF                 = filter(AIF*min_interval,1,Est_IRF_no_Delay,[],2);
% Fit with full 4-parameter 2CXM
conv_result_Larss_with_Delay             = zeros(size(conv_result_no_Delay_IRF));
for i = 1:size(AIF_delay_corrected,1)
    conv_result_Larss_with_Delay(i,:)        = filter(AIF_delay_corrected(i,:)*min_interval,1,fitted_larsson_with_Delay(i,:),[],2);
end
% conv_result_Larss_with_Delay             = filter(AIF_delay_corrected*min_interval,1,fitted_larsson_with_Delay,[],2);
conv_result_Larss_no_Delay               = filter(AIF*min_interval,1,fitted_larsson_no_Delay,[],2);

% Fit with 3-parameter ETM model (special treatment to a delta function)
% conv_result_Larss_with_Delay_High_F      = filter(AIF_delay_corrected*min_interval,1,fitted_larsson_with_Delay_High_F,[],2);
% conv_result_Larss_no_Delay_High_F        = filter(AIF*min_interval,1,fitted_larsson_no_Delay_High_F,[],2);
conv_result_Larss_with_Delay_High_F           = zeros(size(conv_result_no_Delay_IRF));
conv_result_Larss_no_Delay_High_F             = zeros(size(conv_result_no_Delay_IRF));
for i = 1:size(AIF_delay_corrected,1)
    conv_result_Larss_with_Delay_High_F(i,:)       = (Vb_with_Delay_High_F(i) * AIF_delay_corrected(i,:)) + filter(AIF_delay_corrected(i,:)*min_interval,1,fitted_larsson_with_Delay_High_F(i,:),[],2);
    conv_result_Larss_no_Delay_High_F(i,:)         = (Vb_no_Delay_High_F(i)   * AIF                     ) + filter(AIF*min_interval,1,fitted_larsson_no_Delay_High_F(i,:),[],2);
end

% Uptake -> health brain tissues
conv_result_Larss_with_Delay_no_Ve             = zeros(size(conv_result_no_Delay_IRF));
for i = 1:size(AIF_delay_corrected,1)
    conv_result_Larss_with_Delay_no_Ve(i,:)        = filter(AIF_delay_corrected(i,:)*min_interval,1,fitted_larsson_with_Delay_no_Ve(i,:),[],2);
end
% conv_result_Larss_with_Delay_no_Ve       = filter(AIF_delay_corrected*min_interval,1,fitted_larsson_with_Delay_no_Ve,[],2);
conv_result_Larss_no_Delay_no_Ve         = filter(AIF*min_interval,1,fitted_larsson_no_Delay_no_Ve,[],2);

% Uptake, no E -> WM
conv_result_Larss_with_Delay_no_E             = zeros(size(conv_result_no_Delay_IRF));
for i = 1:size(AIF_delay_corrected,1)
    conv_result_Larss_with_Delay_no_E(i,:)        = filter(AIF_delay_corrected(i,:)*min_interval,1,fitted_larsson_with_Delay_no_E(i,:),[],2);
end
% conv_result_Larss_with_Delay_no_E        = filter(AIF_delay_corrected*min_interval,1,fitted_larsson_with_Delay_no_E,[],2);
conv_result_Larss_no_Delay_no_E          = filter(AIF*min_interval,1,fitted_larsson_no_Delay_no_E,[],2);

% Just Vp -> GM  (special treatment to a delta function)
conv_result_Larss_with_Delay_no_E_High_F           = zeros(size(conv_result_no_Delay_IRF));
conv_result_Larss_no_Delay_no_E_High_F             = zeros(size(conv_result_no_Delay_IRF));
for i = 1:size(AIF_delay_corrected,1)
    conv_result_Larss_with_Delay_no_E_High_F(i,:)       = (fitted_larsson_with_Delay_no_E_High_F(1) * AIF_delay_corrected(i,:));
    conv_result_Larss_no_Delay_no_E_High_F(i,:)         = (fitted_larsson_no_Delay_no_E_High_F(1)   * AIF                     );
end
% conv_result_Larss_with_Delay_no_E_High_F = filter(AIF_delay_corrected*min_interval,1,fitted_larsson_with_Delay_no_E_High_F,[],2);
% conv_result_Larss_no_Delay_no_E_High_F   = filter(AIF*min_interval,1,fitted_larsson_no_Delay_no_E_High_F,[],2);
% conv_result_Larss_with_Delay_no_E_High_F = AIF_delay_corrected * fitted_larsson_with_Delay_no_E_High_F(1);
% conv_result_Larss_no_Delay_no_E_High_F   = AIF                 * fitted_larsson_no_Delay_no_E_High_F(1);

conv_result_gaussian                     = filter(AIF*min_interval,1,fitted_gaussian,[],2);
conv_result_double_gaussian              = filter(AIF*min_interval,1,fitted_double_gaussian,[],2);

% Zero negative values
conv_result_no_Delay_IRF(conv_result_no_Delay_IRF<0)                                 = 0;
conv_result_Larss_with_Delay(conv_result_Larss_with_Delay<0)                         = 0;
conv_result_Larss_no_Delay(conv_result_Larss_no_Delay<0)                             = 0;
conv_result_Larss_with_Delay_High_F(conv_result_Larss_with_Delay_High_F<0)           = 0;
conv_result_Larss_no_Delay_High_F(conv_result_Larss_no_Delay_High_F<0)               = 0;
conv_result_Larss_with_Delay_no_Ve(conv_result_Larss_with_Delay_no_Ve<0)             = 0;
conv_result_Larss_no_Delay_no_Ve(conv_result_Larss_no_Delay_no_Ve<0)                 = 0;
conv_result_Larss_with_Delay_no_E(conv_result_Larss_with_Delay_no_E<0)               = 0;
conv_result_Larss_no_Delay_no_E(conv_result_Larss_no_Delay_no_E<0)                   = 0;
conv_result_Larss_with_Delay_no_E_High_F(conv_result_Larss_with_Delay_no_E_High_F<0) = 0;
conv_result_Larss_no_Delay_no_E_High_F(conv_result_Larss_no_Delay_no_E_High_F<0)     = 0;
conv_result_gaussian(conv_result_gaussian<0)                                         = 0;
conv_result_double_gaussian(conv_result_double_gaussian<0)                           = 0;

% Calculate RMS of convolution results comparing to Ct(t)
RMS_ht_no_Delay                  = sqrt( sum( (Ct - conv_result_no_Delay_IRF                ).^2 , 2) );
RMS_Larss_with_Delay             = sqrt( sum( (Ct - conv_result_Larss_with_Delay            ).^2 , 2) );
RMS_Larss_no_Delay               = sqrt( sum( (Ct - conv_result_Larss_no_Delay              ).^2 , 2) );
RMS_Larss_with_Delay_High_F      = sqrt( sum( (Ct - conv_result_Larss_with_Delay_High_F     ).^2 , 2) );
RMS_Larss_no_Delay_High_F        = sqrt( sum( (Ct - conv_result_Larss_no_Delay_High_F       ).^2 , 2) );
RMS_Larss_with_Delay_no_Ve       = sqrt( sum( (Ct - conv_result_Larss_with_Delay_no_Ve      ).^2 , 2) );
RMS_Larss_no_Delay_no_Ve         = sqrt( sum( (Ct - conv_result_Larss_no_Delay_no_Ve        ).^2 , 2) );
RMS_Larss_with_Delay_no_E        = sqrt( sum( (Ct - conv_result_Larss_with_Delay_no_E       ).^2 , 2) );
RMS_Larss_no_Delay_no_E          = sqrt( sum( (Ct - conv_result_Larss_no_Delay_no_E         ).^2 , 2) );
RMS_Larss_with_Delay_no_E_High_F = sqrt( sum( (Ct - conv_result_Larss_with_Delay_no_E_High_F).^2 , 2) );
RMS_Larss_no_Delay_no_E_High_F   = sqrt( sum( (Ct - conv_result_Larss_no_Delay_no_E_High_F  ).^2 , 2) );
RMS_Larss_no_Delay_zero_params   = sqrt( sum( (Ct - zeros(size(Ct))                         ).^2 , 2) );
RMS_gauss                        = sqrt( sum( (Ct - conv_result_gaussian                    ).^2 , 2) );
RMS_double_gauss                 = sqrt( sum( (Ct - conv_result_double_gaussian             ).^2 , 2) );

%% Put all results in return struct
returnStruct = struct();
returnStruct.Flow_with_Delay                          = Flow_with_Delay;
returnStruct.Flow_no_Delay                            = Flow_no_Delay;
returnStruct.Delay_sec_by_Max_Val_with_Delay          = Delay_sec_by_Max_Val_with_Delay;
returnStruct.Delay_sec_by_Max_Val_no_Delay            = Delay_sec_by_Max_Val_no_Delay;
returnStruct.est_delay_by_AIF_correct                 = est_delay_by_AIF_correct;
returnStruct.Est_IRF_with_Delay                       = Est_IRF_with_Delay;
returnStruct.Est_IRF_no_Delay                         = Est_IRF_no_Delay;
returnStruct.fitted_gaussian                          = fitted_gaussian;
returnStruct.conv_result_Larss_with_Delay             = conv_result_Larss_with_Delay;
returnStruct.conv_result_Larss_no_Delay               = conv_result_Larss_no_Delay;
returnStruct.conv_result_Larss_with_Delay_High_F      = conv_result_Larss_with_Delay_High_F;
returnStruct.conv_result_Larss_no_Delay_High_F        = conv_result_Larss_no_Delay_High_F;
returnStruct.conv_result_Larss_with_Delay_no_Ve       = conv_result_Larss_with_Delay_no_Ve;
returnStruct.conv_result_Larss_no_Delay_no_Ve         = conv_result_Larss_no_Delay_no_Ve;
returnStruct.conv_result_Larss_with_Delay_no_E        = conv_result_Larss_with_Delay_no_E;
returnStruct.conv_result_Larss_no_Delay_no_E          = conv_result_Larss_no_Delay_no_E;
returnStruct.conv_result_Larss_with_Delay_no_E_High_F = conv_result_Larss_with_Delay_no_E_High_F;
returnStruct.conv_result_Larss_no_Delay_no_E_High_F   = conv_result_Larss_no_Delay_no_E_High_F;
returnStruct.conv_result_no_Delay_IRF                 = conv_result_no_Delay_IRF;
returnStruct.conv_result_gaussian                     = conv_result_gaussian;
returnStruct.RMS_Larss_with_Delay                     = RMS_Larss_with_Delay;
returnStruct.RMS_Larss_no_Delay                       = RMS_Larss_no_Delay;
returnStruct.RMS_Larss_with_Delay_High_F              = RMS_Larss_with_Delay_High_F;
returnStruct.RMS_Larss_no_Delay_High_F                = RMS_Larss_no_Delay_High_F;
returnStruct.RMS_Larss_with_Delay_no_Ve               = RMS_Larss_with_Delay_no_Ve;
returnStruct.RMS_Larss_no_Delay_no_Ve                 = RMS_Larss_no_Delay_no_Ve;
returnStruct.RMS_Larss_with_Delay_no_E                = RMS_Larss_with_Delay_no_E;
returnStruct.RMS_Larss_no_Delay_no_E                  = RMS_Larss_no_Delay_no_E;
returnStruct.RMS_Larss_with_Delay_no_E_High_F         = RMS_Larss_with_Delay_no_E_High_F;
returnStruct.RMS_Larss_no_Delay_no_E_High_F           = RMS_Larss_no_Delay_no_E_High_F;
returnStruct.RMS_Larss_no_Delay_zero_params           = RMS_Larss_no_Delay_zero_params;
returnStruct.RMS_ht_no_Delay                          = RMS_ht_no_Delay;
returnStruct.RMS_gauss                                = RMS_gauss;
returnStruct.RMS_params_Gauss                         = RMS_params_Gauss;
returnStruct.gaussian_param                           = gaussian_param;
returnStruct.fitted_double_gaussian                   = fitted_double_gaussian;
returnStruct.conv_result_double_gaussian              = conv_result_double_gaussian;
returnStruct.double_gauss_param                       = double_gauss_param;
returnStruct.RMS_double_gauss                         = RMS_double_gauss;
returnStruct.RMS_params_double_gauss                  = RMS_params_double_gauss;
returnStruct.Ktrans_with_Delay                        = Ktrans_with_Delay;
returnStruct.Ktrans_no_Delay                          = Ktrans_no_Delay;
returnStruct.E_with_Delay                             = E_with_Delay;
returnStruct.E_no_Delay                               = E_no_Delay;
returnStruct.E_with_Delay_no_Ve                       = E_with_Delay_no_Ve;
returnStruct.E_no_Delay_no_Ve                         = E_no_Delay_no_Ve;
returnStruct.Ktrans_with_Delay_High_F                 = Ktrans_with_Delay_High_F;
returnStruct.Ktrans_no_Delay_High_F                   = Ktrans_no_Delay_High_F;
returnStruct.Vb_with_Delay                            = Vb_with_Delay;
returnStruct.Vb_no_Delay                              = Vb_no_Delay;
returnStruct.Vb_with_Delay_High_F                     = Vb_with_Delay_High_F;
returnStruct.Vb_no_Delay_High_F                       = Vb_no_Delay_High_F;
returnStruct.Vb_with_Delay_no_Ve                      = Vb_with_Delay_no_Ve;
returnStruct.Vb_no_Delay_no_Ve                        = Vb_no_Delay_no_Ve;
returnStruct.Vb_with_Delay_no_E                       = Vb_with_Delay_no_E;
returnStruct.Vb_no_Delay_no_E                         = Vb_no_Delay_no_E;
returnStruct.Vb_with_Delay_no_E_High_F                = Vb_with_Delay_no_E_High_F;
returnStruct.Vb_no_Delay_no_E_High_F                  = Vb_no_Delay_no_E_High_F;
returnStruct.Ve_with_Delay                            = Ve_with_Delay;
returnStruct.Ve_no_Delay                              = Ve_no_Delay;
returnStruct.Ve_with_Delay_High_F                     = Ve_with_Delay_High_F;
returnStruct.Ve_no_Delay_High_F                       = Ve_no_Delay_High_F;
returnStruct.MTT_with_Delay                           = MTT_with_Delay;
returnStruct.MTT_no_Delay                             = MTT_no_Delay;
returnStruct.Ktrans_Patlak_with_Delay                 = Ktrans_Patlak_with_Delay;
returnStruct.Ktrans_Patlak_no_Delay                   = Ktrans_Patlak_no_Delay;
returnStruct.Vb_Patlak_with_Delay                     = Vb_Patlak_with_Delay;
returnStruct.Vb_Patlak_no_Delay                       = Vb_Patlak_no_Delay;
returnStruct.MTT_Patlak_with_Delay                    = MTT_Patlak_with_Delay;
returnStruct.MTT_Patlak_no_Delay                      = MTT_Patlak_no_Delay;

end

