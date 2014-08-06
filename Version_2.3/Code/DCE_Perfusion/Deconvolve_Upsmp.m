% Upsample data to try and estimate time delay

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

% Create convolution matrix
[ Conv_X_up_samp ]          = Convolution_Matrix( Upsampling_resolution, Sim_AIF_with_noise_T_up_samp );
[ Conv_X_up_samp_no_noise ] = Convolution_Matrix( Upsampling_resolution, Sim_AIF_no_noise_T_up_samp );

% Overwrite B matrix with a new one using the padded time vector
B_mat_upsmp = Create_B_matrix(knots,time_vec_minutes_up_samp,poly_deg-1);

% Deconvolution by regularization for gauss filter
[ridge_regression_gauss_result, b_spline_gauss_result, b_spline_gauss_result_1st_deriv, b_spline_gauss_result_2nd_deriv, b_PCA_gauss_result_2nd_deriv, Sim_Struct.idx_fig]...
    = Regularization_Methods_Simulation(Sim_Ct_gauss_Regul_up_samp,Sim_Ct_gauss_Regul_noise_up_samp,Conv_X_up_samp,Conv_X_up_samp_no_noise,time_vec_minutes_up_samp,...
    lambda_vec_gauss, normalize, Upsampling_resolution, B_mat_upsmp, plot_L_Curve, Sim_Struct.idx_fig , 'Gauss' , Derivative_Time_Devision, plot_flag );

% Deconvolution by regularization for larsson's filter
[ridge_regression_larss_result, b_spline_larss_result, b_spline_larss_result_1st_deriv, b_spline_larss_result_2nd_deriv, b_PCA_gauss_result_2nd_deriv, Sim_Struct.idx_fig]...
    = Regularization_Methods_Simulation(Sim_Ct_larss_Regul_up_samp,Sim_Ct_larss_Regul_noise_up_samp,Conv_X_up_samp,Conv_X_up_samp_no_noise,time_vec_minutes_up_samp,...
    lambda_vec_larss, normalize, Upsampling_resolution, B_mat_upsmp, plot_L_Curve, Sim_Struct.idx_fig , 'Larss' , Derivative_Time_Devision, plot_flag );

