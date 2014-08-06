function [ Return_Struct, idx_fig ] = Estimate_ht_Regularization( Sim_Struct, Verbosity, iter_num, avg_num, idx_fig )

% Take from struct variables used in local function
Sim_Ct_larss_kernel             = Sim_Struct.Sim_Ct_larss_kernel;
Sim_Ct_gauss_kernel             = Sim_Struct.Sim_Ct_gauss_kernel_noise;
min_interval                    = Sim_Struct.min_interval;
time_vec_minutes                = Sim_Struct.time_vec_minutes;
gauss_filter                    = Sim_Struct.gauss_filter;
Sim_AIF_no_noise                = Sim_Struct.Sim_AIF_no_noise;
Sim_AIF_with_noise              = Sim_Struct.Sim_AIF_with_noise;
Sim_Ct_gauss_kernel_noise       = Sim_Struct.Sim_Ct_gauss_kernel_noise;
Sim_Ct_larss_kernel_noise       = Sim_Struct.Sim_Ct_larss_kernel_noise;
Use_Cyclic_Conv_4_ht_est        = Sim_Struct.Use_Cyclic_Conv_4_ht_est;
Use_Upsampling_Delay_Comp       = Sim_Struct.Use_Upsampling_Delay_Comp;
lambda_vec_gauss                = Sim_Struct.lambda_vec_gauss;
normalize                       = Sim_Struct.normalize;
B_mat                           = Sim_Struct.B_mat;
PCA_B_mat                       = Sim_Struct.PCA_B_mat;
plot_L_Curve                    = Sim_Struct.plot_L_Curve;
Derivative_Time_Devision        = Sim_Struct.Derivative_Time_Devision;
plot_flag                       = Sim_Struct.plot_flag;
lambda_vec_larss                = Sim_Struct.lambda_vec_larss;

if strcmp(Verbosity,'Full')
    display('-I- Estimating h(t) using Regularization...');
end

Sim_AIF_no_noise_Regul   = Sim_AIF_no_noise(:,iter_num,avg_num);
Sim_AIF_with_noise_Regul = Sim_AIF_with_noise(:,iter_num,avg_num);
Sim_Ct_gauss_Regul       = Sim_Ct_gauss_kernel(:,iter_num,avg_num);
Sim_Ct_gauss_Regul_noise = Sim_Ct_gauss_kernel_noise(:,iter_num,avg_num);
Sim_Ct_larss_Regul       = Sim_Ct_larss_kernel(:,iter_num,avg_num);
Sim_Ct_larss_Regul_noise = Sim_Ct_larss_kernel_noise(:,iter_num,avg_num);

% Create convolution indices
[ Conv_X ]          = Convolution_Matrix( min_interval(iter_num), Sim_AIF_with_noise_Regul );
[ Conv_X_no_noise ] = Convolution_Matrix( min_interval(iter_num), Sim_AIF_no_noise_Regul );

Sim_Struct.Sim_AIF_no_noise_Regul   = Sim_AIF_no_noise_Regul;
Sim_Struct.Sim_AIF_with_noise_Regul = Sim_AIF_with_noise_Regul;

% Trying to emphasize the first knots
%knots   = [time_vec_minutes(1:4) time_vec_minutes(5:4:end)  ];

%end_idx = round(max(size(time_vec_minutes))/5);
%knots   = time_vec_minutes(1:1:end_idx);
%knots   = time_vec_minutes(1:1:11);

if Use_Cyclic_Conv_4_ht_est
    % Cyclic de-convolution
    [ridge_regression_larss_result, b_spline_larss_result, b_spline_larss_result_1st_deriv, b_spline_larss_result_2nd_deriv, b_PCA_larss_result_2nd_deriv,...
     ridge_regression_gauss_result, b_spline_gauss_result, b_spline_gauss_result_1st_deriv, b_spline_gauss_result_2nd_deriv, b_PCA_gauss_result_2nd_deriv,...
     idx_fig] = Cyclic_Deconvolve( Sim_Struct, Verbosity, iter_num, avg_num, idx_fig );
    
elseif Use_Upsampling_Delay_Comp
    % Up-sample, no cyclic de-convolution
    Deconvolve_Upsmp;
    
else
    % Regular de-convolution, no correction for delay
    
    % Deconvolution by regularization for gauss filter
    [ridge_regression_gauss_result, b_spline_gauss_result, b_spline_gauss_result_1st_deriv, b_spline_gauss_result_2nd_deriv, b_PCA_gauss_result_2nd_deriv, idx_fig]...
        = Regularization_Methods_Simulation(Sim_Ct_gauss_Regul, Sim_Ct_gauss_Regul_noise,Conv_X,Conv_X_no_noise,time_vec_minutes,...
        lambda_vec_gauss, normalize, min_interval(iter_num), B_mat, PCA_B_mat, plot_L_Curve, idx_fig , 'Gauss' , Derivative_Time_Devision, plot_flag );
    
    % Deconvolution by regularization for larsson's filter
    [ridge_regression_larss_result, b_spline_larss_result, b_spline_larss_result_1st_deriv, b_spline_larss_result_2nd_deriv, b_PCA_larss_result_2nd_deriv, idx_fig]...
        = Regularization_Methods_Simulation(Sim_Ct_larss_Regul, Sim_Ct_larss_Regul_noise,Conv_X,Conv_X_no_noise,time_vec_minutes,...
        lambda_vec_larss, normalize, min_interval(iter_num), B_mat, PCA_B_mat, plot_L_Curve, idx_fig , 'Larss' , Derivative_Time_Devision, plot_flag );
    
end

% Return vars
Return_Struct.Sim_AIF_no_noise_Regul           = Sim_AIF_no_noise_Regul;
Return_Struct.Sim_AIF_with_noise_Regul         = Sim_AIF_with_noise_Regul;
Return_Struct.ridge_regression_gauss_result    = ridge_regression_gauss_result;
Return_Struct.ridge_regression_larss_result    = ridge_regression_larss_result;
Return_Struct.b_spline_gauss_result            = b_spline_gauss_result;
Return_Struct.b_spline_gauss_result_1st_deriv  = b_spline_gauss_result_1st_deriv;
Return_Struct.b_spline_gauss_result_2nd_deriv  = b_spline_gauss_result_2nd_deriv;
Return_Struct.b_PCA_gauss_result_2nd_deriv     = b_PCA_gauss_result_2nd_deriv;
Return_Struct.Sim_Ct_larss_Regul               = Sim_Ct_larss_Regul;
Return_Struct.Sim_Ct_larss_Regul_noise         = Sim_Ct_larss_Regul_noise;
Return_Struct.Conv_X_no_noise                  = Conv_X_no_noise;
Return_Struct.idx_fig                          = idx_fig;

% The following can be modified in AIF_Delay_Correct
Return_Struct.b_spline_larss_result            = b_spline_larss_result;
Return_Struct.b_spline_larss_result_1st_deriv  = b_spline_larss_result_1st_deriv;
Return_Struct.b_spline_larss_result_2nd_deriv  = b_spline_larss_result_2nd_deriv;
Return_Struct.b_PCA_larss_result_2nd_deriv     = b_PCA_larss_result_2nd_deriv;

end