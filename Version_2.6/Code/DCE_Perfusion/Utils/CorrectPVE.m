function [ AIF_estimated_ICA_Corrected, Scale_Factor ] = CorrectPVE( AIF_estimated_ICA, Vein_estimated_ICA, Sim_Struct )

% Take needed parameters from struct
min_interval                = Sim_Struct.min_interval;
lambda_vec_larss            = Sim_Struct.lambda_vec_larss;
normalize                   = Sim_Struct.normalize;
plot_L_Curve                = Sim_Struct.plot_L_Curve;
idx_fig                     = Sim_Struct.idx_fig;
Derivative_Time_Devision    = Sim_Struct.Derivative_Time_Devision;
plot_flag                   = Sim_Struct.plot_flag;
RealData_Flag               = Sim_Struct.RealData_Flag;
time_vec_minutes            = Sim_Struct.time_vec_minutes;
poly_deg                    = Sim_Struct.poly_deg;
knots                       = Sim_Struct.knots;

% % Regular integral method
% Integral_Art                = trapz(time_vec_minutes,AIF_estimated_ICA);
% Integral_Vein               = trapz(time_vec_minutes,Vein_estimated_ICA);
% Scale_Factor                = Integral_Vein / Integral_Art;
% AIF_estimated_ICA_Corrected = Scale_Factor * AIF_estimated_ICA;

% Deconvolution method
[ Conv_Matrix ]             = Convolution_Matrix( min_interval, AIF_estimated_ICA );
Conv_Matrix_no_noise        = Conv_Matrix;
filter_type                 = 'Larss';
B_PCA                       = zeros(size(Conv_Matrix));
B_mat                       = Create_B_matrix(knots, time_vec_minutes, poly_deg-1);

[ ~, ~, ~, b_spline_result_2nd_deriv, ~, ~, ~, ~ ] =  ...
    Regularization_Methods_Simulation( Vein_estimated_ICA', AIF_estimated_ICA', Conv_Matrix, Conv_Matrix_no_noise, time_vec_minutes, lambda_vec_larss, normalize, min_interval, B_mat, B_PCA, plot_L_Curve, idx_fig, filter_type, Derivative_Time_Devision, plot_flag, RealData_Flag );

Scale_Factor                = 1./ trapz(time_vec_minutes,b_spline_result_2nd_deriv);
AIF_estimated_ICA_Corrected = Scale_Factor * AIF_estimated_ICA;


end



