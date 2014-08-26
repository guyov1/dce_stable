function [ Final_Filter_Estimation_Larss ] = Choose_Needed_Ht( Filter_Est_Chosen, est_larss_filter_Wiener_noise, ridge_regression_larss_result, b_spline_larss_result, b_spline_larss_result_1st_deriv, b_spline_larss_result_2nd_deriv, b_PCA_larss_result, b_PCA_larss_result_1st_deriv, b_PCA_larss_result_2nd_deriv)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Choose which filter estimation to use
switch Filter_Est_Chosen
    case 'Wiener'
        Final_Filter_Estimation_Larss = est_larss_filter_Wiener_noise';
    case 'Ridge'
        Final_Filter_Estimation_Larss = ridge_regression_larss_result;
    case 'Spline'
        Final_Filter_Estimation_Larss = b_spline_larss_result;
    case 'Spline_1st'
        Final_Filter_Estimation_Larss = b_spline_larss_result_1st_deriv;
    case 'Spline_2nd'
        Final_Filter_Estimation_Larss = b_spline_larss_result_2nd_deriv;
    case 'PCA'
        Final_Filter_Estimation_Larss = b_PCA_larss_result;
    case 'PCA_1st'
        Final_Filter_Estimation_Larss = b_PCA_larss_result_1st_deriv;
    case 'PCA_2nd'
        Final_Filter_Estimation_Larss = b_PCA_larss_result_2nd_deriv;
    otherwise
        error('Unrecognized filter estimation method!');
end

end

