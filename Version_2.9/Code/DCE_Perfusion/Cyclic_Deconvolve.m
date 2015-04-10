function [ridge_regression_larss_result, b_spline_larss_result, b_spline_larss_result_1st_deriv, b_spline_larss_result_2nd_deriv, b_PCA_larss_result, b_PCA_larss_result_1st_deriv, b_PCA_larss_result_2nd_deriv,ridge_regression_gauss_result, b_spline_gauss_result, b_spline_gauss_result_1st_deriv, b_spline_gauss_result_2nd_deriv, b_PCA_gauss_result, b_PCA_gauss_result_1st_deriv, b_PCA_gauss_result_2nd_deriv, idx_fig] = Cyclic_Deconvolve( Sim_Struct, Verbosity, iter_num, avg_num, idx_fig )

% Take from struct variables used in local function
Sim_Ct_larss_kernel              = Sim_Struct.Sim_Ct_larss_kernel(:,iter_num,avg_num);
Sim_Ct_gauss_Regul               = Sim_Struct.Sim_Ct_gauss_kernel(:,iter_num,avg_num);
Sim_Ct_gauss_Regul_noise         = Sim_Struct.Sim_Ct_gauss_kernel(:,iter_num,avg_num);
Sim_Ct_larss_Regul               = Sim_Struct.Sim_Ct_larss_kernel(:,iter_num,avg_num);
Sim_Ct_larss_Regul_noise         = Sim_Struct.Sim_Ct_larss_kernel(:,iter_num,avg_num);
min_interval                     = Sim_Struct.min_interval;
gauss_filter                     = Sim_Struct.gauss_filter(:,iter_num);
larss_filter                     = Sim_Struct.larss_filter(:,iter_num);
Sim_AIF_with_noise               = Sim_Struct.Sim_AIF_with_noise(:,iter_num,avg_num);
Sim_Ct_larss_kernel_noise        = Sim_Struct.Sim_Ct_larss_kernel_noise(:,iter_num,avg_num);
lambda_vec_gauss                 = Sim_Struct.lambda_vec_gauss;
normalize                        = Sim_Struct.normalize;
plot_L_Curve                     = Sim_Struct.plot_L_Curve;
Derivative_Time_Devision         = Sim_Struct.Derivative_Time_Devision;
plot_flag                        = Sim_Struct.plot_flag;
lambda_vec_larss                 = Sim_Struct.lambda_vec_larss;
Sim_AIF_delayed_no_noise         = Sim_Struct.Sim_AIF_delayed_no_noise(:,iter_num,avg_num);
Sim_AIF_no_noise_Regul           = Sim_Struct.Sim_AIF_no_noise_Regul;
Sim_AIF_with_noise_Regul         = Sim_Struct.Sim_AIF_with_noise_Regul;
total_sim_time_min               = Sim_Struct.total_sim_time_min;
knots                            = Sim_Struct.knots;
poly_deg                         = Sim_Struct.poly_deg;
Use_Upsampling_and_Cyclic        = Sim_Struct.Use_Upsampling_and_Cyclic;
RealData_Flag                    = Sim_Struct.RealData_Flag;

% Cyclic de-convolution

% Outputs the deconvolved h(t), using cylic de-convolution and an optional
% upsampling

%% ---------------- Try to apply cyclic convolution to compensate for delay --------------

% Set padding length (total length = total length + Pad_times* ( total length )
AIF_len            = length(Sim_AIF_with_noise_Regul);

if Sim_Struct.Cyclic_End_Padding
    Sim_Ct_gauss_Regul_padded       = [Sim_Ct_gauss_Regul ;zeros(AIF_len + 1,1)];
    Sim_Ct_gauss_Regul_noise_padded = [Sim_Ct_gauss_Regul_noise ;zeros(AIF_len + 1,1)];
    Sim_Ct_larss_Regul_padded       = [Sim_Ct_larss_Regul ;zeros(AIF_len + 1,1)];
    Sim_Ct_larss_Regul_noise_padded = [Sim_Ct_larss_Regul_noise ;zeros(AIF_len + 1,1)];
else
    Sim_Ct_gauss_Regul_padded       = [zeros(AIF_len + 1,1) ; Sim_Ct_gauss_Regul ];
    Sim_Ct_gauss_Regul_noise_padded = [zeros(AIF_len + 1,1) ; Sim_Ct_gauss_Regul_noise];
    Sim_Ct_larss_Regul_padded       = [zeros(AIF_len + 1,1) ; Sim_Ct_larss_Regul ];
    Sim_Ct_larss_Regul_noise_padded = [zeros(AIF_len + 1,1) ; Sim_Ct_larss_Regul_noise];
end

% Time vector for AIF and Ct(t)
time_vec_minutes_padded    = (0:2*AIF_len).* min_interval(iter_num);

% Create circular convolution matrix
[ Conv_X_cyclic ]          = Cyclic_Convolution_Matrix( min_interval(iter_num), Sim_AIF_with_noise_Regul);
[ Conv_X_cyclic_no_noise ] = Cyclic_Convolution_Matrix( min_interval(iter_num), Sim_AIF_no_noise_Regul);

% Overwrite B matrix with a new one using the padded time vector
B_mat_cyclic          = Create_B_matrix(knots,time_vec_minutes_padded,poly_deg-1);
B_PCA_cyclic          = PCA_basis(Sim_Struct, time_vec_minutes_padded);
% Take # of eigen-vectors similar to B-splines
num_cols_B_mat_cyclic = size(B_mat_cyclic,2);
B_PCA_cyclic          = B_PCA_cyclic(:,1:num_cols_B_mat_cyclic);

% Deconvolution by regularization for gauss filter
[ridge_regression_gauss_result, b_spline_gauss_result, b_spline_gauss_result_1st_deriv, b_spline_gauss_result_2nd_deriv, b_PCA_gauss_result, b_PCA_gauss_result_1st_deriv, b_PCA_gauss_result_2nd_deriv, Sim_Struct.idx_fig]...
    = Regularization_Methods_Simulation(Sim_Ct_gauss_Regul_padded,Sim_Ct_gauss_Regul_noise_padded,Conv_X_cyclic,Conv_X_cyclic_no_noise,time_vec_minutes_padded,...
    lambda_vec_gauss, normalize, min_interval(iter_num), B_mat_cyclic, B_PCA_cyclic, plot_L_Curve, Sim_Struct.idx_fig , 'Gauss' , Derivative_Time_Devision, plot_flag, RealData_Flag );

% Deconvolution by regularization for larsson's filter
[ridge_regression_larss_result, b_spline_larss_result, b_spline_larss_result_1st_deriv, b_spline_larss_result_2nd_deriv, b_PCA_larss_result, b_PCA_larss_result_1st_deriv, b_PCA_larss_result_2nd_deriv, Sim_Struct.idx_fig]...
    = Regularization_Methods_Simulation(Sim_Ct_larss_Regul_padded,Sim_Ct_larss_Regul_noise_padded,Conv_X_cyclic,Conv_X_cyclic_no_noise,time_vec_minutes_padded,...
    lambda_vec_larss, normalize, min_interval(iter_num), B_mat_cyclic, B_PCA_cyclic,  plot_L_Curve, Sim_Struct.idx_fig , 'Larss' , Derivative_Time_Devision, plot_flag, RealData_Flag );

% Remove the padding out of the h_t result
if Sim_Struct.Cyclic_End_Padding
    ridge_regression_gauss_result   = ridge_regression_gauss_result(1:AIF_len);
    b_spline_gauss_result           = b_spline_gauss_result(1:AIF_len);
    b_spline_gauss_result_1st_deriv = b_spline_gauss_result_1st_deriv(1:AIF_len);
    b_spline_gauss_result_2nd_deriv = b_spline_gauss_result_2nd_deriv(1:AIF_len);
    b_PCA_gauss_result              = b_PCA_gauss_result(1:AIF_len);
    b_PCA_gauss_result_1st_deriv    = b_PCA_gauss_result_1st_deriv(1:AIF_len);
    b_PCA_gauss_result_2nd_deriv    = b_PCA_gauss_result_2nd_deriv(1:AIF_len);
    
    ridge_regression_larss_result   = ridge_regression_larss_result(1:AIF_len);
    b_spline_larss_result           = b_spline_larss_result(1:AIF_len);
    b_spline_larss_result_1st_deriv = b_spline_larss_result_1st_deriv(1:AIF_len);
    b_spline_larss_result_2nd_deriv = b_spline_larss_result_2nd_deriv(1:AIF_len);
    b_PCA_larss_result              = b_PCA_larss_result(1:AIF_len);
    b_PCA_larss_result_1st_deriv    = b_PCA_larss_result_1st_deriv(1:AIF_len);
    b_PCA_larss_result_2nd_deriv    = b_PCA_larss_result_2nd_deriv(1:AIF_len);
else
    ridge_regression_gauss_result   = ridge_regression_gauss_result(AIF_len+1:end);
    b_spline_gauss_result           = b_spline_gauss_result(AIF_len+1:end);
    b_spline_gauss_result_1st_deriv = b_spline_gauss_result_1st_deriv(AIF_len+1:end);
    b_spline_gauss_result_2nd_deriv = b_spline_gauss_result_2nd_deriv(AIF_len+1:end);
    b_PCA_gauss_result              = b_PCA_gauss_result(AIF_len+1:end);
    b_PCA_gauss_result_1st_deriv    = b_PCA_gauss_result_1st_deriv(AIF_len+1:end);
    b_PCA_gauss_result_2nd_deriv    = b_PCA_gauss_result_2nd_deriv(AIF_len+1:end);
    
    ridge_regression_larss_result   = ridge_regression_larss_result(AIF_len+1:end);
    b_spline_larss_result           = b_spline_larss_result(AIF_len+1:end);
    b_spline_larss_result_1st_deriv = b_spline_larss_result_1st_deriv(AIF_len+1:end);
    b_spline_larss_result_2nd_deriv = b_spline_larss_result_2nd_deriv(AIF_len+1:end);
    b_PCA_larss_result              = b_PCA_larss_result(AIF_len+1:end);
    b_PCA_larss_result_1st_deriv    = b_PCA_larss_result_1st_deriv(AIF_len+1:end);
    b_PCA_larss_result_2nd_deriv    = b_PCA_larss_result_2nd_deriv(AIF_len+1:end);
end

%% -------------------------------------

if (Use_Upsampling_and_Cyclic)
    
    % Set upsampling factor according to paramaters
    UpSampFactor                     = round(min_interval / Upsampling_resolution) ;
    % Time vector for AIF and Ct(t)
    time_vec_minutes_up_samp         = (0:(num_time_stamps*UpSampFactor) - 1).* Upsampling_resolution;
    
    % ------------------   TO REMOVE!!! -----------------------------------
    time_vec_minutes                       = time_vec_minutes_up_samp;
    Est_Wiener_Filter_real_PSD_larss       = interp(Est_Wiener_Filter_real_PSD_larss,UpSampFactor);
    larss_filter                           = interp(larss_filter,UpSampFactor);
    Est_Wiener_Filter_est_PSD_larss        = interp(Est_Wiener_Filter_est_PSD_larss,UpSampFactor);
    Sim_AIF_with_noise                     = interp(Sim_AIF_with_noise,UpSampFactor);
    Sim_Ct_larss_kernel                    = interp(Sim_Ct_larss_kernel,UpSampFactor);
    Sim_AIF_delayed_no_noise               = interp(Sim_AIF_delayed_no_noise,UpSampFactor);
    est_larss_filter_Wiener_noise          = interp(est_larss_filter_Wiener_noise,UpSampFactor);
    Sim_Ct_larss_kernel_noise              = interp(Sim_Ct_larss_kernel_noise,UpSampFactor);
    % ---------------------------------------------------------------------
    
    % Upsample Cts and AIF
    Sim_Ct_gauss_Regul_up_samp       = interp(Sim_Ct_gauss_Regul,UpSampFactor);
    Sim_Ct_gauss_Regul_noise_up_samp = interp(Sim_Ct_gauss_Regul_noise,UpSampFactor);
    Sim_Ct_larss_Regul_up_samp       = interp(Sim_Ct_larss_Regul,UpSampFactor);
    Sim_Ct_larss_Regul_noise_up_samp = interp(Sim_Ct_larss_Regul_noise,UpSampFactor);
    Sim_AIF_with_noise_T_up_samp     = interp(Sim_AIF_with_noise_Regul,UpSampFactor);
    Sim_AIF_no_noise_T_up_samp       = interp(Sim_AIF_no_noise_Regul,UpSampFactor);
    
    % Set padding length (total length = total length + Pad_times* ( total length )
    
    Pad_times = 1;
    AIF_upsmp_len                     = length(Sim_AIF_with_noise_T_up_samp);
    
    Sim_Ct_gauss_Regul_up_samp        = [Sim_Ct_gauss_Regul_up_samp ;zeros(AIF_upsmp_len + 1,1)];
    Sim_Ct_gauss_Regul_noise_up_samp  = [Sim_Ct_gauss_Regul_noise_up_samp ;zeros(AIF_upsmp_len + 1,1)];
    Sim_Ct_larss_Regul_up_samp        = [Sim_Ct_larss_Regul_up_samp ;zeros(AIF_upsmp_len + 1,1)];
    Sim_Ct_larss_Regul_noise_up_samp  = [Sim_Ct_larss_Regul_noise_up_samp ;zeros(AIF_upsmp_len + 1,1)];
    
    % Time vector for AIF and Ct(t)
    time_vec_minutes_up_samp           = (0:2*AIF_upsmp_len).* Upsampling_resolution;
    
    % Create circular convolution matrix
    [ Conv_X_up_samp ]          = Cyclic_Convolution_Matrix( Upsampling_resolution, Sim_AIF_with_noise_T_up_samp);
    [ Conv_X_up_samp_no_noise ] = Cyclic_Convolution_Matrix( Upsampling_resolution, Sim_AIF_no_noise_T_up_samp);
    
    % Overwrite B matrix with a new one using the padded time vector
    B_mat_cyclic_upsmp    = Create_B_matrix(knots,time_vec_minutes_up_samp,poly_deg-1);
    B_PCA_cyclic_upsmp    = PCA_basis(Sim_Struct, time_vec_minutes_up_samp);
    % Take # of eigen-vectors similar to B-splines
    num_cols_B_mat_cyclic = size(B_mat_cyclic_upsmp,2);
    B_PCA_cyclic_upsmp    = B_PCA_cyclic_upsmp(:,1:num_cols_B_mat_cyclic);
    
    % Deconvolution by regularization for gauss filter
    [ridge_regression_gauss_result, b_spline_gauss_result, b_spline_gauss_result_1st_deriv, b_spline_gauss_result_2nd_deriv, b_PCA_gauss_result, b_PCA_gauss_result_1st_deriv, b_PCA_gauss_result_2nd_deriv, Sim_Struct.idx_fig]...
        = Regularization_Methods_Simulation(Sim_Ct_gauss_Regul_up_samp,Sim_Ct_gauss_Regul_noise_up_samp,Conv_X_up_samp,Conv_X_up_samp_no_noise,time_vec_minutes_up_samp,...
        lambda_vec_gauss, normalize, Upsampling_resolution, B_mat_cyclic_upsmp, B_PCA_cyclic_upsmp, plot_L_Curve, Sim_Struct.idx_fig , 'Gauss' , Derivative_Time_Devision, plot_flag, RealData_Flag );
    
    % Deconvolution by regularization for larsson's filter
    [ridge_regression_larss_result, b_spline_larss_result, b_spline_larss_result_1st_deriv, b_spline_larss_result_2nd_deriv, b_PCA_larss_result, b_PCA_larss_result_1st_deriv, b_PCA_larss_result_2nd_deriv, Sim_Struct.idx_fig]...
        = Regularization_Methods_Simulation(Sim_Ct_larss_Regul_up_samp,Sim_Ct_larss_Regul_noise_up_samp,Conv_X_up_samp,Conv_X_up_samp_no_noise,time_vec_minutes_up_samp,...
        lambda_vec_larss, normalize, Upsampling_resolution, B_mat_cyclic_upsmp, B_PCA_cyclic_upsmp, plot_L_Curve, Sim_Struct.idx_fig , 'Larss' , Derivative_Time_Devision, plot_flag, RealData_Flag );
    
    % Remove the padding out of the h_t result
    
    b_PCA_larss_result_2nd_deriv    = b_PCA_larss_result_2nd_deriv(1:AIF_upsmp_len);
    
    % Remove the padding out of the h_t result
    if Sim_Struct.Cyclic_End_Padding
        ridge_regression_gauss_result   = ridge_regression_gauss_result(1:AIF_upsmp_len);
        b_spline_gauss_result           = b_spline_gauss_result(1:AIF_upsmp_len);
        b_spline_gauss_result_1st_deriv = b_spline_gauss_result_1st_deriv(1:AIF_upsmp_len);
        b_spline_gauss_result_2nd_deriv = b_spline_gauss_result_2nd_deriv(1:AIF_upsmp_len);
        b_PCA_gauss_result              = b_PCA_gauss_result(1:AIF_upsmp_len);
        b_PCA_gauss_result_1st_deriv    = b_PCA_gauss_result_1st_deriv(1:AIF_upsmp_len);
        b_PCA_gauss_result_2nd_deriv    = b_PCA_gauss_result_2nd_deriv(1:AIF_upsmp_len);
        
        ridge_regression_larss_result   = ridge_regression_larss_result(1:AIF_upsmp_len);
        b_spline_larss_result           = b_spline_larss_result(1:AIF_upsmp_len);
        b_spline_larss_result_1st_deriv = b_spline_larss_result_1st_deriv(1:AIF_upsmp_len);
        b_spline_larss_result_2nd_deriv = b_spline_larss_result_2nd_deriv(1:AIF_upsmp_len);
        b_PCA_larss_result              = b_PCA_larss_result(1:AIF_upsmp_len);
        b_PCA_larss_result_1st_deriv    = b_PCA_larss_result_1st_deriv(1:AIF_upsmp_len);
        
    else
        ridge_regression_gauss_result   = ridge_regression_gauss_result(AIF_upsmp_len+1:end);
        b_spline_gauss_result           = b_spline_gauss_result(AIF_upsmp_len+1:end);
        b_spline_gauss_result_1st_deriv = b_spline_gauss_result_1st_deriv(AIF_upsmp_len+1:end);
        b_spline_gauss_result_2nd_deriv = b_spline_gauss_result_2nd_deriv(AIF_upsmp_len+1:end);
        b_PCA_gauss_result              = b_PCA_gauss_result(AIF_upsmp_len+1:end);
        b_PCA_gauss_result_1st_deriv    = b_PCA_gauss_result_1st_deriv(AIF_upsmp_len+1:end);
        b_PCA_gauss_result_2nd_deriv    = b_PCA_gauss_result_2nd_deriv(AIF_upsmp_len+1:end);
        
        ridge_regression_larss_result   = ridge_regression_larss_result(AIF_upsmp_len+1:end);
        b_spline_larss_result           = b_spline_larss_result(AIF_upsmp_len+1:end);
        b_spline_larss_result_1st_deriv = b_spline_larss_result_1st_deriv(AIF_upsmp_len+1:end);
        b_spline_larss_result_2nd_deriv = b_spline_larss_result_2nd_deriv(AIF_upsmp_len+1:end);
        b_PCA_larss_result              = b_PCA_larss_result(AIF_upsmp_len+1:end);
        b_PCA_larss_result_1st_deriv    = b_PCA_larss_result_1st_deriv(AIF_upsmp_len+1:end);
        b_PCA_larss_result_2nd_deriv    = b_PCA_larss_result_2nd_deriv(AIF_upsmp_len+1:end);
    end
    
    
    
end

end