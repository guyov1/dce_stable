function [ Flow_vec_with_Delay, Flow_vec_no_Delay, Delay_sec_by_Max_Val_with_Delay, Delay_sec_by_Max_Val_no_Delay, est_delay_by_AIF_correct, t_delay_single_gauss_sec,...
    sigma_seconds_single_gauss, Amp_single_gauss, Est_IRF_with_Delay, Est_IRF_no_Delay, fitted_gaussian, ...
    conv_result_Larss_with_Delay, conv_result_Larss_no_Delay, conv_result_Larss_no_Delay_High_F, conv_result_Larss_no_Delay_no_E, conv_result_Larss_no_Delay_no_E_High_F, ...
    conv_result_no_Delay_IRF, conv_result_gaussian, ...
    RMS_Larss_with_Delay, RMS_Larss_no_Delay, RMS_Larss_with_Delay_High_F, RMS_Larss_no_Delay_High_F, RMS_Larss_with_Delay_no_E_3D, ...
    RMS_Larss_no_Delay_no_E, RMS_Larss_with_Delay_no_E_High_F, RMS_Larss_no_Delay_no_E_High_F, RMS_Larss_no_Delay_zero_params,...
    RMS_ht_no_Delay, RMS_gauss, RMS_params_Gauss, fitted_double_gaussian, conv_result_double_gaussian, ...
    double_gaussian_param_vec, RMS_double_gauss, RMS_params_double_gauss, Ktrans_with_Delay_vec, Ktrans_no_Delay_vec, ...
    E_with_Delay, E_no_Delay, Ktrans_with_Delay_High_F, Ktrans_no_Delay_High_F, ...
    Vb_with_Delay, Vb_no_Delay, Vb_with_Delay_High_F, Vb_no_Delay_High_F, Vb_with_Delay_no_E, Vb_no_Delay_no_E, Vb_with_Delay_no_E_High_F, Vb_no_Delay_no_E_High_F, ...
    Ve_with_Delay, Ve_no_Delay, Ve_with_Delay_High_F, Ve_no_Delay_High_F, ...
    MTT_with_Delay, MTT_no_Delay, Ktrans_Patlak_with_Delay_vec, Ktrans_Patlak_no_Delay_vec, ...
    Vb_Patlak_with_Delay_vec, Vb_Patlak_no_Delay_vec, MTT_Patlak_with_Delay_vec, MTT_Patlak_no_Delay_vec ] ...
    = Get_Ht_Deconvolving(Sim_Struct, AIF, Ct , Output_directory, Subject_name, Force_RealData_Calc, Verbosity)

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
%   - RMS_params_Gauss, the root mean square error of the gaussian using the parameters
%          and h(t)  (if too big, h(t) is not really a gaussian)
%   - calculated_double_gaussian, the fitted double gaussian to ht
%   - conv_result_double_gaussian,  the result of AIF convolved with double gaussian
%   - double_gaussian_param_vec, all the double gaussian parameters (means, vars and amplitudes).
%   - RMS_double_gauss, the root mean square error of conv_result with estimated double gaussian and Ct.
%   - RMS_params_double_gauss, the root mean square error of the double gaussian using the parameters
%                   and h(t)  (if too big, h(t) is not really double gaussian)
%   - Ktrans   - estimation for permeability according to Larsson
%   - Vb   - estimation for blood volume according to Larsson
%   - Ve   - estimation for EES volume according to Larsson
%   - MTT  - estimation for Mean Transient Time according to Larsson
%   - Ktrans_Patlak_vec  - estimation for permeability according to Patlak
%   - Vb_Patlak_vec  - estimation for blood volume according to Patlak
%   - MTT_Patlak_vec - estimation for Mean Transient Time according to Patlak


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
    
    Sim_Ct_T                 = NaN;
    %Conv_Matrix              = ; % Already Assigned
    Conv_Matrix_no_noise     = NaN;
    %time_vec_minutes         = ; % Already Assigned
    %lambda_vec               = ; % Already Assigned
    %normalize                = ; % Already Assigned
    %min_interval             = ; % Already Assigned
    %B_mat                    = ; % Already Assigned
    B_PCA                    = NaN;
    %plot_L_Curve             = ; % Already Assigned
    idx_fig                  = 1;
    filter_type              = 'Larss';
    %Derivative_Time_Devision = ; % Already Assigned
    plot_flag                = false;
    RealData_Flag            = Sim_Struct.RealData_Flag;
    
    % --------------- AIF delay correction parameters ------------------
    In_Struct1                                    = struct;
    In_Struct1.normalize                          = normalize;
    In_Struct1.B_mat                              = B_mat;
    In_Struct1.PCA_B_mat                          = NaN;
    In_Struct1.plot_L_Curve                       = plot_L_Curve;
    In_Struct1.Derivative_Time_Devision           = Derivative_Time_Devision;
    In_Struct1.lambda_vec_larss                   = lambda_vec_larss;
    In_Struct1.min_interval                       = min_interval;
    In_Struct1.time_vec_minutes                   = time_vec_minutes;
    In_Struct1.Upsampling_resolution              = Upsampling_resolution;
    In_Struct1.Max_Time_Delay                     = Max_Time_Delay;
    In_Struct1.Min_Time_Delay                     = Min_Time_Delay;
    In_Struct1.Use_Upsampling_Delay_Comp          = Use_Upsampling_Delay_Comp;
    In_Struct1.LowerBound_Larsson                 = LowerBound_Larsson;
    In_Struct1.UpperBound_Larsson                 = UpperBound_Larsson;
    In_Struct1.algorithm_options                  = algorithm_options;
    In_Struct1.Hct                                = Hct; % Try to read it from patient
    In_Struct1.RMS_Smooth_Around_Bolus            = RMS_Smooth_Around_Bolus;
    In_Struct1.RMS_Smooth                         = RMS_Smooth;
    In_Struct1.Diff_From_Bolus                    = Diff_From_Bolus;
    In_Struct1.additional_AIF_delay_sec           = 0;
    In_Struct1.BiExp2CTC_RMS_Ratio                = BiExp2CTC_RMS_Ratio;
    In_Struct1.plot_flag                          = false;
    In_Struct1.Adjusted_Larsson_Model             = Adjusted_Larsson_Model;
    In_Struct1.Filter_Est_Chosen                  = Filter_Est_Chosen;
    In_Struct1.Vb_low                             = Vb_low;
    In_Struct1.RealData_Flag                      = RealData_Flag;
    In_Struct1.Simple_AIF_Delay_Correct           = Simple_AIF_Delay_Correct;
    In_Struct1.Patlak_Est_Type                    = Patlak_Est_Type;
    In_Struct1.Ktrans                             = NaN;       % Simulation ground truth values
    In_Struct1.Vb_larss                           = NaN; % Simulation ground truth values
    In_Struct1.init_Ve_guess                      = Sim_Struct.init_Ve_guess;
    In_Struct1.FMS_Algorithm                      = Sim_Struct.FMS_Algorithm;
    In_Struct1.LQ_Model_AIF_Delay_Correct         = Sim_Struct.LQ_Model_AIF_Delay_Correct;
    In_Struct1.filter_type                        = Sim_Struct.filter_type;

    In_Struct2                                    = struct;
    In_Struct2.Sim_AIF_with_noise_Regul           = AIF;
    In_Struct2.Sim_Ct_larss_Regul                 = NaN;
    In_Struct2.Conv_X_no_noise                    = NaN;
    iter_num                                      = 1;
    avg_num                                       = 1;
    % ------------------------------------------------------------
    
    
    if(exist(Mat_File_Ht,'file') && ~Force_RealData_Calc)
        load(Mat_File_Ht);
        display('--------------------------------------------------------');
        display('-I- Starting h(t) estimation using regularization...');
        display('--------------------------------------------------------');
    else
        
        
        display('--------------------------------------------------------');
        display('-I- Starting h(t) estimation using regularization...');
        display('--------------------------------------------------------');
        
        AIF_delay_corrected = zeros(size(Ct));
        
        %for j=1:num_voxels
        parfor j=1:num_voxels
            tic;
            
            [ ~, ~, ~, b_spline_result_2nd_deriv_no_Delay, ~, ~, ~, ~ ] =  ...
                Regularization_Methods_Simulation( Sim_Ct_T, Ct(j,:)', Conv_Matrix, Conv_Matrix_no_noise, time_vec_minutes, lambda_vec_larss, normalize, min_interval, B_mat, B_PCA, plot_L_Curve, idx_fig, filter_type, Derivative_Time_Devision, plot_flag, RealData_Flag );
            
            Est_IRF_no_Delay(j,:) = b_spline_result_2nd_deriv_no_Delay;
            
            % Estimate delay
            % Correct h(t) estimation if it seems we have delay in AIF
            if Correct_estimation_due_to_delay
                [est_delay_by_AIF_correct(j), AIF_delay_corrected(j ,:), Est_IRF_with_Delay(j,:), ~] = ...
                    AIF_Delay_Correct(In_Struct1, In_Struct2, Est_IRF_no_Delay(j,:), Ct(j,:)', Verbosity, iter_num, avg_num, idx_fig);
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
        save(Mat_File_Ht, 'Est_IRF_with_Delay','Est_IRF_no_Delay', 'est_delay_by_AIF_correct', 'AIF_delay_corrected');
        
    end
    
end

%% Estimate parameters for all voxels by curve fitting


% Check if already estimated parameters for all voxels
Mat_File_Perfusion_Parameters = [Output_directory 'Estimated_Perfusion_Parameters_' Subject_name '.mat'];

if(exist(Mat_File_Perfusion_Parameters,'file') && ~Force_RealData_Calc)
    load(Mat_File_Perfusion_Parameters);
else
    
    if Parallel_Real_Data_Est
        [ ...
          Flow_vec_with_Delay, Flow_vec_no_Delay, Delay_sec_by_Max_Val_with_Delay, Delay_sec_by_Max_Val_no_Delay, ...
          t_delay_single_gauss_sec, sigma_seconds_single_gauss, Amp_single_gauss, ...
          fitted_larsson_with_Delay, fitted_larsson_no_Delay, fitted_larsson_with_Delay_High_F, fitted_larsson_no_Delay_High_F, ...
          fitted_larsson_with_Delay_no_E, fitted_larsson_no_Delay_no_E, fitted_larsson_with_Delay_no_E_High_F, fitted_larsson_no_Delay_no_E_High_F, ...
          fitted_gaussian, fitted_double_gaussian, double_gaussian_param_vec, ...
          Ktrans_with_Delay_vec, Ktrans_no_Delay_vec, ...
          E_with_Delay, E_no_Delay, Ktrans_with_Delay_High_F, Ktrans_no_Delay_High_F,  ...
          Vb_with_Delay, Vb_no_Delay, Vb_with_Delay_High_F, Vb_no_Delay_High_F, ...
          Vb_with_Delay_no_E, Vb_no_Delay_no_E, Vb_with_Delay_no_E_High_F, Vb_no_Delay_no_E_High_F, ...
          Ve_with_Delay, Ve_no_Delay, Ve_with_Delay_High_F, Ve_no_Delay_High_F, MTT_with_Delay, MTT_no_Delay, ...
          Ktrans_Patlak_with_Delay_vec, Ktrans_Patlak_no_Delay_vec, Vb_Patlak_with_Delay_vec, Vb_Patlak_no_Delay_vec, ...
          MTT_Patlak_with_Delay_vec, MTT_Patlak_no_Delay_vec ] ...
          = Parallel_Params_Est_Real_Data(Sim_Struct, Est_IRF_with_Delay, Est_IRF_no_Delay, Ct, AIF_delay_corrected,AIF, idx_fig );
    else
        
        [ ...
          Flow_vec_with_Delay, Flow_vec_no_Delay, Delay_sec_by_Max_Val_with_Delay, Delay_sec_by_Max_Val_no_Delay, ...
          t_delay_single_gauss_sec, sigma_seconds_single_gauss, Amp_single_gauss, ...
          fitted_larsson_with_Delay, fitted_larsson_no_Delay, fitted_larsson_with_Delay_High_F, fitted_larsson_no_Delay_High_F, ...
          fitted_larsson_with_Delay_no_E, fitted_larsson_no_Delay_no_E, fitted_larsson_with_Delay_no_E_High_F, fitted_larsson_no_Delay_no_E_High_F, ...
          fitted_gaussian, fitted_double_gaussian, double_gaussian_param_vec, ...
          Ktrans_with_Delay_vec, Ktrans_no_Delay_vec, ...
          E_with_Delay, E_no_Delay, Ktrans_with_Delay_High_F, Ktrans_no_Delay_High_F,  ...
          Vb_with_Delay, Vb_no_Delay, Vb_with_Delay_High_F, Vb_no_Delay_High_F, ...
          Vb_with_Delay_no_E, Vb_no_Delay_no_E, Vb_with_Delay_no_E_High_F, Vb_no_Delay_no_E_High_F, ...
          Ve_with_Delay, Ve_no_Delay, Ve_with_Delay_High_F, Ve_no_Delay_High_F, MTT_with_Delay, MTT_no_Delay, ...
          Ktrans_Patlak_with_Delay_vec, Ktrans_Patlak_no_Delay_vec, Vb_Patlak_with_Delay_vec, Vb_Patlak_no_Delay_vec, ...
          MTT_Patlak_with_Delay_vec, MTT_Patlak_no_Delay_vec ] ...
          = Serial_Params_Est_Real_Data(Sim_Struct, Est_IRF_with_Delay, Est_IRF_no_Delay, Ct, AIF_delay_corrected, AIF, idx_fig );
    end
    
    save(Mat_File_Perfusion_Parameters,'Flow_vec_with_Delay','Flow_vec_no_Delay','t_delay_single_gauss_sec','sigma_seconds_single_gauss',...
        'Amp_single_gauss', 'Ktrans_with_Delay_vec','Ktrans_no_Delay_vec', 'E_with_Delay', 'E_no_Delay', 'Ktrans_with_Delay_High_F', 'Ktrans_no_Delay_High_F', ...
        'Vb_with_Delay', 'Vb_no_Delay', 'Vb_with_Delay_High_F', 'Vb_no_Delay_High_F', 'Vb_with_Delay_no_E', 'Vb_no_Delay_no_E', 'Vb_with_Delay_no_E_High_F', 'Vb_no_Delay_no_E_High_F', ...
        'Ve_with_Delay', 'Ve_no_Delay', 'Ve_with_Delay_High_F', 'Ve_no_Delay_High_F', 'MTT_with_Delay', 'MTT_no_Delay',...
        'Ktrans_Patlak_with_Delay_vec','Ktrans_Patlak_no_Delay_vec', 'Vb_Patlak_with_Delay_vec','Vb_Patlak_no_Delay_vec',...
        'MTT_Patlak_with_Delay_vec','MTT_Patlak_no_Delay_vec','Delay_sec_by_Max_Val_with_Delay',...
        'Delay_sec_by_Max_Val_no_Delay','double_gaussian_param_vec','fitted_gaussian','fitted_double_gaussian', ...
        'fitted_larsson_with_Delay','fitted_larsson_no_Delay', 'fitted_larsson_with_Delay_High_F', 'fitted_larsson_no_Delay_High_F', 'fitted_larsson_no_Delay_no_E', 'fitted_larsson_with_Delay_no_E', 'fitted_larsson_with_Delay_no_E_High_F', 'fitted_larsson_no_Delay_no_E_High_F');
    
end

% Calculate RMS of convolution result comparing to Ct(t)
RMS_params_Gauss        = sqrt( sum( (fitted_gaussian        - Est_IRF_no_Delay).^2 ,2) );
RMS_params_double_gauss = sqrt( sum( (fitted_double_gaussian - Est_IRF_no_Delay).^2 ,2) );

%% Filter AIF through kernel

% Filter the AIF with the estimated ht
%conv_result_ht       = filter(Est_ht,1,AIF);
if Sim_Struct.ignore_time_delta
    conv_result_no_Delay_IRF                 = filter(AIF,1,Est_IRF_no_Delay,[],2);
    conv_result_Larss_with_Delay             = filter(AIF,1,fitted_larsson_with_Delay,[],2);
    conv_result_Larss_no_Delay               = filter(AIF,1,fitted_larsson_no_Delay,[],2);
    conv_result_Larss_with_Delay_High_F      = filter(AIF,1,fitted_larsson_with_Delay_High_F,[],2);
    conv_result_Larss_no_Delay_High_F        = filter(AIF,1,fitted_larsson_no_Delay_High_F,[],2);
    conv_result_Larss_with_Delay_no_E        = filter(AIF,1,fitted_larsson_with_Delay_no_E,[],2);
    conv_result_Larss_no_Delay_no_E          = filter(AIF,1,fitted_larsson_no_Delay_no_E,[],2);
    conv_result_Larss_with_Delay_no_E_High_F = filter(AIF,1,fitted_larsson_with_Delay_no_E_High_F,[],2); 
    conv_result_Larss_no_Delay_no_E_High_F   = filter(AIF,1,fitted_larsson_no_Delay_no_E_High_F,[],2);
else
    conv_result_no_Delay_IRF                 = filter(AIF*min_interval,1,Est_IRF_no_Delay,[],2);
    conv_result_Larss_with_Delay             = filter(AIF*min_interval,1,fitted_larsson_with_Delay,[],2);
    conv_result_Larss_no_Delay               = filter(AIF*min_interval,1,fitted_larsson_no_Delay,[],2);
    conv_result_Larss_with_Delay_High_F      = filter(AIF*min_interval,1,fitted_larsson_with_Delay_High_F,[],2);
    conv_result_Larss_no_Delay_High_F        = filter(AIF*min_interval,1,fitted_larsson_no_Delay_High_F,[],2);
    conv_result_Larss_with_Delay_no_E        = filter(AIF*min_interval,1,fitted_larsson_with_Delay_no_E,[],2);
    conv_result_Larss_no_Delay_no_E          = filter(AIF*min_interval,1,fitted_larsson_no_Delay_no_E,[],2);
    conv_result_Larss_with_Delay_no_E_High_F = filter(AIF*min_interval,1,fitted_larsson_with_Delay_no_E_High_F,[],2); 
    conv_result_Larss_no_Delay_no_E_High_F   = filter(AIF*min_interval,1,fitted_larsson_no_Delay_no_E_High_F,[],2);
end

% Filter the AIF with the gaussian kernel
%conv_result_gaussian = filter(calculated_gaussian*min_interval,1,AIF);
if Sim_Struct.ignore_time_delta
    conv_result_gaussian        = filter(AIF,1,fitted_gaussian,[],2);
    conv_result_double_gaussian = filter(AIF,1,fitted_double_gaussian,[],2);
else
    conv_result_gaussian        = filter(AIF*min_interval,1,fitted_gaussian,[],2);
    conv_result_double_gaussian = filter(AIF*min_interval,1,fitted_double_gaussian,[],2);
end

% Zero negative values
conv_result_no_Delay_IRF(conv_result_no_Delay_IRF<0)                                 = 0;
conv_result_Larss_with_Delay(conv_result_Larss_with_Delay<0)                         = 0;
conv_result_Larss_no_Delay(conv_result_Larss_no_Delay<0)                             = 0;
conv_result_Larss_with_Delay_High_F(conv_result_Larss_with_Delay_High_F<0)           = 0;
conv_result_Larss_no_Delay_High_F(conv_result_Larss_no_Delay_High_F<0)               = 0;
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
RMS_Larss_with_Delay_no_E_3D     = sqrt( sum( (Ct - conv_result_Larss_with_Delay_no_E       ).^2 , 2) );
RMS_Larss_no_Delay_no_E          = sqrt( sum( (Ct - conv_result_Larss_no_Delay_no_E         ).^2 , 2) );
RMS_Larss_with_Delay_no_E_High_F = sqrt( sum( (Ct - conv_result_Larss_with_Delay_no_E_High_F).^2 , 2) );
RMS_Larss_no_Delay_no_E_High_F   = sqrt( sum( (Ct - conv_result_Larss_no_Delay_no_E_High_F  ).^2 , 2) );
%RMS_Larss_no_Delay_zero_params = sqrt( sum( (Ct - repmat(mean(Ct,2), 1, size(Ct,2))).^2 , 2) );
RMS_Larss_no_Delay_zero_params   = sqrt( sum( (Ct - zeros(size(Ct))                         ).^2 , 2) );

RMS_gauss                        = sqrt( sum( (Ct - conv_result_gaussian).^2, 2) );
RMS_double_gauss                 = sqrt( sum( (Ct - conv_result_double_gaussian).^2, 2) );

% Change all places of with 0 to Inf (because 0 is unrealistic)

% RMS_ht_no_Delay                  (RMS_ht_no_Delay                  == 0) = Inf;
% RMS_Larss_with_Delay             (RMS_Larss_with_Delay             == 0) = Inf;
% RMS_Larss_no_Delay               (RMS_Larss_no_Delay               == 0) = Inf;
% RMS_Larss_with_Delay_High_F      (RMS_Larss_with_Delay_High_F      == 0) = Inf;
% RMS_Larss_no_Delay_High_F        (RMS_Larss_no_Delay_High_F        == 0) = Inf;
% RMS_Larss_with_Delay_no_E_3D     (RMS_Larss_with_Delay_no_E_3D     == 0) = Inf;
% RMS_Larss_no_Delay_no_E          (RMS_Larss_no_Delay_no_E          == 0) = Inf;
% RMS_Larss_with_Delay_no_E_High_F (RMS_Larss_with_Delay_no_E_High_F == 0) = Inf;
% RMS_Larss_no_Delay_no_E_High_F   (RMS_Larss_no_Delay_no_E_High_F   == 0) = Inf;
% %RMS_Larss_no_Delay_zero_params = sqrt( sum( (Ct - repmat(mean(Ct,2), 1, size(Ct,2))).^2 , 2) );
% RMS_Larss_no_Delay_zero_params   (RMS_Larss_no_Delay_zero_params   == 0) = Inf;
% 
% RMS_gauss                        (RMS_gauss                        == 0) = Inf;
% RMS_double_gauss                 (RMS_double_gauss                 == 0) = Inf;

end

