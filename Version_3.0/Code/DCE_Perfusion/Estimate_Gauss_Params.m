function [est_sigma_noise, est_t_d_noise, est_amp_noise,  idx_fig ] = Estimate_Gauss_Params( Sim_Struct, ht_Struct , Verbosity, iter_num, avg_num, idx_fig)

% Take from struct variables used in local function
min_interval                    = Sim_Struct.min_interval;
time_vec_minutes                = Sim_Struct.time_vec_minutes;
gauss_filter                    = Sim_Struct.gauss_filter;
Sim_AIF_with_noise              = Sim_Struct.Sim_AIF_with_noise;
Sim_Ct_gauss_kernel_noise       = Sim_Struct.Sim_Ct_gauss_kernel_noise;
plot_flag                       = Sim_Struct.plot_flag;
Filter_Est_Chosen               = Sim_Struct.Filter_Est_Chosen;
FMS_Algorithm                   = Sim_Struct.FMS_Algorithm;
Init_Guess_Gaussian             = Sim_Struct.Init_Guess_Gaussian;
LowerBound_Gauss                = Sim_Struct.LowerBound_Gauss;
UpperBound_Gauss                = Sim_Struct.UpperBound_Gauss;
algorithm_options               = Sim_Struct.algorithm_options;
est_gauss_filter_Wiener_noise   = Sim_Struct.est_gauss_filter_Wiener_noise;
ridge_regression_gauss_result   = ht_Struct.ridge_regression_gauss_result;
b_spline_gauss_result           = ht_Struct.b_spline_gauss_result;
b_spline_gauss_result_1st_deriv = ht_Struct.b_spline_gauss_result_1st_deriv;
b_spline_gauss_result_2nd_deriv = ht_Struct.b_spline_gauss_result_2nd_deriv;
b_PCA_gauss_result_1st_deriv    = ht_Struct.b_PCA_gauss_result_1st_deriv;
b_PCA_gauss_result_2nd_deriv    = ht_Struct.b_PCA_gauss_result_2nd_deriv;

Final_Filter_Estimation_Gauss   = ht_struct.Final_Filter_Estimation_Gauss;
  

if strcmp(Verbosity,'Full')
    display('-I- Estimating Gaussian parameters...');
end

% lsqcurvefit parameters are:
%      analytic function, initial parameters, time vector, data points ,lower
%      and upper bounds and algorithm options

% The analytic funcational of a gaussian function
Gaussian_function = @(x,t) Gaussian( t,x(1),x(2),x(3) );

if strcmp(FMS_Algorithm,'trust-region-reflective')
    [est_params_Gauss_noise,residue_norm_Gauss_noise,residual_Gauss_noise,exitflag_Gauss_noise,algo_info_Gauss_noise] = ...
        lsqcurvefit(Gaussian_function,Init_Guess_Gaussian,time_vec_minutes',Final_Filter_Estimation_Gauss,LowerBound_Gauss,UpperBound_Gauss,algorithm_options);
    
elseif strcmp(FMS_Algorithm,'levenberg-marquardt')
    [est_params_Gauss_noise,residue_norm_Gauss_noise,residual_Gauss_noise,exitflag_Gauss_noise,algo_info_Gauss_noise] = ...
        lsqcurvefit(Gaussian_function,Init_Guess_Gaussian,time_vec_minutes',Final_Filter_Estimation_Gauss,[],[],algorithm_options);
    
end

% Assign variables with Gaussian parameters estimation
est_t_d_noise   = est_params_Gauss_noise(1);
est_var_noise   = est_params_Gauss_noise(2);
est_amp_noise   = est_params_Gauss_noise(3);
est_sigma_noise = sqrt(est_var_noise);

% Plot results
if (plot_flag)
    
    fig_num = figure;
    subplot(3,1,1);
    hold on;
    h1  = plot(time_vec_minutes,gauss_filter,'g');
    h2  = plot(time_vec_minutes,gauss_filter,'g*');
    h3  = plot(time_vec_minutes,est_gauss_filter_Wiener_noise,'b');
    h4  = plot(time_vec_minutes,est_gauss_filter_Wiener_noise,'b+');
    h5  = plot(time_vec_minutes,ridge_regression_gauss_result,'k');
    h6  = plot(time_vec_minutes,ridge_regression_gauss_result,'ko');
    h7  = plot(time_vec_minutes,b_spline_gauss_result,'c');
    h8  = plot(time_vec_minutes,b_spline_gauss_result,'cx');
    h9  = plot(time_vec_minutes,b_spline_gauss_result_1st_deriv,'m');
    h10 = plot(time_vec_minutes,b_spline_gauss_result_1st_deriv,'md');
    h11 = plot(time_vec_minutes,b_spline_gauss_result_2nd_deriv,'r');
    h12 = plot(time_vec_minutes,b_spline_gauss_result_2nd_deriv,'rs');
    h13 = plot(time_vec_minutes,b_PCA_larss_result_2nd_deriv ,'y');
    h14 = plot(time_vec_minutes,b_PCA_larss_result_2nd_deriv ,'yh');
    
    title('True and estimated h(t) - Gauss','FontWeight','bold');
    xlabel('Time [Min]');
    hold off;
    legend([h2 h4 h6 h8 h10 h12 h14],'Orig. h(t)',...
        'Est. h(t) - Wiener','Est. h(t) - Ridge','Est. h(t) - Spline','Est. h(t) - Spline 1st deriv','Est. h(t) - Spline 2nd deriv','Est. h(t) - PCA 2nd deriv');
    
    subplot(3,1,2);
    hold on;
    h1  = plot(time_vec_minutes,Sim_Ct_gauss_kernel_noise,'g');
    h2  = plot(time_vec_minutes,Sim_Ct_gauss_kernel_noise,'g*');
    h3  = plot(time_vec_minutes,filter(est_gauss_filter_Wiener_noise*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num)),'b');
    h4  = plot(time_vec_minutes,filter(est_gauss_filter_Wiener_noise*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num)),'b+');
    h5  = plot(time_vec_minutes,filter(ridge_regression_gauss_result*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num)),'k');
    h6  = plot(time_vec_minutes,filter(ridge_regression_gauss_result*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num)),'ko');
    h7  = plot(time_vec_minutes,filter(b_spline_gauss_result*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num)),'c');
    h8  = plot(time_vec_minutes,filter(b_spline_gauss_result*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num)),'cx');
    h9  = plot(time_vec_minutes,filter(b_spline_gauss_result_1st_deriv*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num)),'m');
    h10 = plot(time_vec_minutes,filter(b_spline_gauss_result_1st_deriv*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num)),'md');
    h11 = plot(time_vec_minutes,filter(b_spline_gauss_result_2nd_deriv*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num)),'r');
    h12 = plot(time_vec_minutes,filter(b_spline_gauss_result_2nd_deriv*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num)),'rs');
    h13 = plot(time_vec_minutes,filter(b_PCA_larss_result_2nd_deriv*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num)),'y');
    h14 = plot(time_vec_minutes,filter(b_PCA_larss_result_2nd_deriv*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num)),'yh');
    
    title('True and estimated Ct(t)','FontWeight','bold');
    xlabel('Time [Min]');
    hold off;
    legend([h2 h4 h6 h8 h10 h12 h14],'Orig. Ct(t)','Est. Ct(t) - Wiener',...
        'Est. Ct(t) - Ridge','Est. Ct(t) - Spline','Est. Ct(t) - Spline 1st deriv','Est. Ct(t) - Spline 2nd deriv','Est. Ct(t) - PCA 2nd deriv');
    
    subplot(3,1,3);
    filter_from_estimated_params = Gaussian( time_vec_minutes, est_t_d_noise, est_var_noise, est_amp_noise );
    plot(time_vec_minutes,gauss_filter,'b*',time_vec_minutes,filter_from_estimated_params,'r.');
    %title('h(t)- Estimated params (Blue).Orig h(t) (Red dots)');
    title('Estimated Params fit result with noise','FontWeight','bold');
    legend('True h(t)','h(t) - Est. Params');
    xlabel('Time [Min]');
    
    % Print result to PDF
    [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Filters_Gauss.png', './Run_Output/',...
        'Estimated Filters after de-convolution - Gaussian', 'EstFiltersGauss');
    
end

end