function [ Return_Struct, idx_fig ] = Estimate_Larss_Params( Sim_Struct, ht_Struct, est_t_d_noise, est_larss_filter_Wiener_noise, Verbosity, iter_num, avg_num, idx_fig)

% Take from struct variables used in local function
Filter_Est_Chosen               = Sim_Struct.Filter_Est_Chosen;
FMS_Algorithm                   = Sim_Struct.FMS_Algorithm;
Init_Guess_Gaussian             = Sim_Struct.Init_Guess_Gaussian;
LowerBound_Gauss                = Sim_Struct.LowerBound_Gauss;
UpperBound_Gauss                = Sim_Struct.UpperBound_Gauss;
Sim_AIF_delayed_no_noise        = Sim_Struct.Sim_AIF_delayed_no_noise;
Ki                              = Sim_Struct.Ki;
F                               = Sim_Struct.F;
E                               = Sim_Struct.E;
PS                              = Sim_Struct.PS;
Vd                              = Sim_Struct.Vd;
Vb_larss                        = Sim_Struct.Vb_larss;
Ve_larss                        = Sim_Struct.Ve_larss;
MTT                             = Sim_Struct.MTT;
SNR_single                      = Sim_Struct.SNR_single;
Ignore_Gaussian_Calculation     = Sim_Struct.Ignore_Gaussian_Calculation;
Sim_Ct_larss_kernel             = Sim_Struct.Sim_Ct_larss_kernel;
Sim_AIF_with_noise              = Sim_Struct.Sim_AIF_with_noise;
Sim_Ct_larss_kernel_noise       = Sim_Struct.Sim_Ct_larss_kernel_noise;
plot_flag                       = Sim_Struct.plot_flag;
min_interval                    = Sim_Struct.min_interval(iter_num);
time_vec_minutes                = Sim_Struct.time_vec_minutes;
LowerBound_Larsson              = Sim_Struct.LowerBound_Larsson;
UpperBound_Larsson              = Sim_Struct.UpperBound_Larsson;
algorithm_options               = Sim_Struct.algorithm_options;
Hct                             = Sim_Struct.Hct;
additional_AIF_delay_sec        = Sim_Struct.additional_AIF_delay_sec;
larss_filter                    = Sim_Struct.larss_filter;
Correct_estimation_due_to_delay = Sim_Struct.Correct_estimation_due_to_delay;
est_delay_by_AIF_correct        = Sim_Struct.est_delay_by_AIF_correct;
Check_Sourbron_Estimate         = Sim_Struct.Check_Sourbron_Estimate;
Random_init_F_guess             = Sim_Struct.Random_init_F_guess;

ridge_regression_larss_result   = ht_Struct.ridge_regression_larss_result;
b_spline_larss_result           = ht_Struct.b_spline_larss_result;
b_spline_larss_result_1st_deriv = ht_Struct.b_spline_larss_result_1st_deriv;
b_spline_larss_result_2nd_deriv = ht_Struct.b_spline_larss_result_2nd_deriv;
b_PCA_larss_result_2nd_deriv    = ht_Struct.b_PCA_larss_result_2nd_deriv;

   
    
if strcmp(Verbosity,'Full')
    display('-I- Estimating Larsson parameters...');
end

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
        Final_Filter_Estimation_Larss = b_PCA_larss_result_2nd_deriv;
    otherwise
        error('Unrecognized filter estimation method!');
end


% F estimated is the maximum value of F*IRF
est_F_noise         = max(Final_Filter_Estimation_Larss);

% Delay of the AIF will be calculated according to the place of the
% maximum value of F*IRF
max_index           = (Final_Filter_Estimation_Larss == est_F_noise);

% Translate from minutes
if (Correct_estimation_due_to_delay)
    est_Delay_sec_noise =  est_delay_by_AIF_correct;
else
    est_Delay_sec_noise =  time_vec_minutes(max_index) * 60;
end


%% Patalk Estimation
[est_Ki_Patlak_noise, est_Vb_Patlak_noise ,E_Patlak_est, idx_fig] = Patlak_Estimation(Sim_Struct, est_F_noise, Verbosity, iter_num, avg_num, idx_fig);

%% Using Murase estimation for Vb,Ki - Alternative to Patlak (Using estimated F)

% -------------- Vb estimation ---------------
%Sim_Ct_larss_kernel = filter(larss_filter*min_interval,1,Sim_AIF_delayed_no_noise);
Y_vec_Vb            = ( est_F_noise* cumtrapz(time_vec_minutes,Sim_Ct_larss_kernel) ) ./ ...
    ( est_F_noise* cumtrapz(time_vec_minutes,Sim_AIF_delayed_no_noise) - Sim_Ct_larss_kernel);

% Take samples after bolus to better estimate Vb
[~, bolus_idx]      = max(diff(Sim_AIF_with_noise(:,iter_num,avg_num)));
diff_from_bolus_min = 40/60; % The difference in seconds(minutes) from the bolus to look on
diff_from_bolus_idx = round(diff_from_bolus_min/min_interval);
stable_idx          = (bolus_idx + diff_from_bolus_idx) : length(Sim_AIF_with_noise(:,iter_num,avg_num));
Vb_murase_estimate  = mean(Y_vec_Vb(stable_idx));

% -------------- Ki estimation ---------------
Y_vec_Ki            = ( Vb_murase_estimate*est_F_noise*cumtrapz(time_vec_minutes,Sim_AIF_delayed_no_noise) - ...
    - Sim_Ct_larss_kernel - est_F_noise* cumtrapz(time_vec_minutes,Sim_Ct_larss_kernel) ) ./ ...
    ( cumtrapz(time_vec_minutes,Sim_Ct_larss_kernel));

Ki_murase_estimate  = mean(Y_vec_Ki(stable_idx));


if (plot_flag)
    
    fig_num = figure;
    
    subplot(2,1,1);
    hold on;
    h1 = plot(time_vec_minutes,Y_vec_Vb,'k*');
    h2 = plot(time_vec_minutes(stable_idx),Y_vec_Vb(stable_idx),'og');
    hold off;
    title('Vb Murase Estimation','FontWeight','bold');
    xlabel('Time [Min]');
    legend([h1 h2],'Vb vector','Stable points for estimation');
    
    %             % Display Larsson's parameters
    %             annotation('textbox',...
    %                 [0.7 0.17 0.3 0.1],...
    %                 'String',{['Vb = ' num2str(Vb_larss) '  [mL/100g]'],...
    %                 ['Vb Murase = ' num2str(Vb_murase_estimate,'%.2f') '  [mL/100g]']},...
    %                 'FontSize',8,...
    %                 'FontName','Arial',...
    %                 'LineStyle','-',...
    %                 'EdgeColor',[0 0 0],...
    %                 'LineWidth',2,...
    %                 'BackgroundColor',[1 1 1],...
    %                 'Color',[0.0 0.0 0]);
    
    subplot(2,1,2);
    hold on;
    h1 = plot(time_vec_minutes,Y_vec_Ki,'k*');
    h2 = plot(time_vec_minutes(stable_idx),Y_vec_Ki(stable_idx),'og');
    hold off;
    title('Ki Murase Estimation','FontWeight','bold');
    xlabel('Time [Min]');
    legend([h1 h2],'Ki vector','Stable points for estimation');
    
    % Display Larsson's parameters
    annotation('textbox',...
        [0.7 0.09 0.3 0.1],...
        'String',{['Ki = ' num2str(Ki) '  [mL/100g]'],...
        ['Ki Murase = ' num2str(Ki_murase_estimate,'%.2f') '  [mL/100g]'],...
        ['Vb = ' num2str(Vb_larss) '  [mL/100g]'],...
        ['Vb Murase = ' num2str(Vb_murase_estimate,'%.2f') '  [mL/100g]']},...
        'FontSize',8,...
        'FontName','Arial',...
        'LineStyle','-',...
        'EdgeColor',[0 0 0],...
        'LineWidth',2,...
        'BackgroundColor',[1 1 1],...
        'Color',[0.0 0.0 0]);
    
    % Print result to PDF
    [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Vb_Ki_Murase_Estimate.png', './Run_Output/',...
        'Murase Estimate for Vb, Ki', 'VbKiMuraseEst');
    
end

%% Estimate using two compartment model using initial Patlak's guess

% Initial Guess for non-linear curve fitting for Larsson (Vb, E, Ve)
Init_Guess_Larsson = [ est_Vb_Patlak_noise E_Patlak_est 5];
%Init_Guess_Larsson = [50 0.15 10];

% The analytic funcational of a Larsson function
Larsson_function             = @(x,t) Larsson_Filter( t, est_F_noise, x(1), x(2), x(3), Hct(iter_num));
Larsson_function_4_Sourbron  = @(x,t) Larsson_Filter( t, x(1)       , x(2), x(3), x(4), Hct(iter_num));

% lsqcurvefit parameters are:
% analytic function, initial parameters, time vector, data points ,lower
% and upper bounds and algorithm options


if strcmp(FMS_Algorithm,'trust-region-reflective')
    
    [est_params_Larsson_noise,residue_norm_Larsson_noise,residual_Larsson_noise,exitflag_Larsson_noise,algo_info_Larsson_noise] = ...
        lsqcurvefit(Larsson_function,Init_Guess_Larsson,time_vec_minutes',Final_Filter_Estimation_Larss/est_F_noise,...
        LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
    
    % Estimate gaussian parameters when Larss filter is used
    if ~Ignore_Gaussian_Calculation
        
        [est_params_Larss_Usiung_Gauss_noise,~,~,~,~] = ...
            lsqcurvefit(Gaussian_function,Init_Guess_Gaussian,time_vec_minutes',Final_Filter_Estimation_Larss,LowerBound_Gauss,UpperBound_Gauss,algorithm_options);
    end
    
    if Check_Sourbron_Estimate
        
        % Gues initial F randomly
        if (Random_init_F_guess)
            F_init_guess        = rand*100;
        else
            F_init_guess        = est_F_noise;
        end
        
        Init_Guess_Sourbron = [F_init_guess Init_Guess_Larsson];
        
        F_low               = 0;
        F_max               = 150; 
        LowerBound_Sourbron = [F_low LowerBound_Larsson];
        UpperBound_Sourbron = [F_max UpperBound_Larsson];
        
        % F, Vb, E, Ve
        [est_params_Sourbron_noise,~,~,~,~] = ...
        lsqcurvefit(Larsson_function_4_Sourbron,Init_Guess_Sourbron,time_vec_minutes',Final_Filter_Estimation_Larss/F_init_guess,...
        LowerBound_Sourbron,UpperBound_Sourbron,algorithm_options);
        
        
    end
    
elseif strcmp(FMS_Algorithm,'levenberg-marquardt')
    
    [est_params_Larsson_noise,residue_norm_Larsson_noise,residual_Larsson_noise,exitflag_Larsson_noise,algo_info_Larsson_noise] = ...
        lsqcurvefit(Larsson_function,Init_Guess_Larsson,time_vec_minutes',Final_Filter_Estimation_Larss/est_F_noise,...
        [],[],algorithm_options);
    
    % Estimate gaussian parameters when Larss filter is used
    if ~Ignore_Gaussian_Calculation
        
        [est_params_Larss_Usiung_Gauss_noise,~,~,~,~] = ...
            lsqcurvefit(Gaussian_function,Init_Guess_Gaussian,time_vec_minutes',Final_Filter_Estimation_Larss,[],[],algorithm_options);
        
    end
    
    if Check_Sourbron_Estimate
        
        % Gues initial F randomly
        if (Random_init_F_guess)
            F_init_guess        = rand*100;
        else
            F_init_guess        = est_F_noise;
        end
        
        Init_Guess_Sourbron = [F_init_guess Init_Guess_Larsson];
        
        % F, Vb, E, Ve
        [est_params_Sourbron_noise,~,~,~,~] = ...
        lsqcurvefit(Larsson_function,Init_Guess_Sourbron,time_vec_minutes',Final_Filter_Estimation_Larss/F_init_guess,...
        [],[],algorithm_options);
        
    end
    
end

% Estimate AIF delay using Larss Filter with Gauss de-convolution
if ~Ignore_Gaussian_Calculation
    est_Delay_sec_using_Gaussian_noise =  est_params_Larss_Usiung_Gauss_noise(1) * 60;
else
    est_Delay_sec_using_Gaussian_noise = NaN;
end

% Assigning two compartment parameters estimation
est_Vb_Two_Comp_noise           = est_params_Larsson_noise(1);
E_Two_Comp_est                  = est_params_Larsson_noise(2);
est_Ki_Two_Comp_noise           = E_Two_Comp_est * est_F_noise;
Ve_Two_Comp_est                 = est_params_Larsson_noise(3);

% Assigning two compartment parameters estimation using Sourbron
if Check_Sourbron_Estimate
    est_F_Two_Comp_Sourbron_noise               = est_params_Sourbron_noise(1);
    est_Vb_Two_Comp_Sourbron_noise              = est_params_Sourbron_noise(2);
    E_Two_Comp_Sourbron_est                     = est_params_Sourbron_noise(3);
    est_Ki_Two_Comp_Sourbron_noise              = E_Two_Comp_Sourbron_est * est_F_Two_Comp_Sourbron_noise;
    Ve_Two_Comp_Sourbron_est                    = est_params_Sourbron_noise(4);
end

% Estimate MTT (also in normal tissue)
est_IRF                         = Final_Filter_Estimation_Larss / est_F_noise;
est_MTT_noise                   = cumtrapz(time_vec_minutes,est_IRF);
est_MTT_noise                   = est_MTT_noise(end);
est_MTT_normal_tis_noise        =  est_Vb_Patlak_noise / est_F_noise;
% Estimate Vd (also in normal tissue)
est_Vd_noise                    = est_F_noise * est_MTT_noise;
est_Vd_normal_tis_noise         = est_F_noise * est_MTT_normal_tis_noise;
% Estimate E
est_E_noise                     = E_Two_Comp_est;
% Estimate PS
est_PS_noise                    = -est_F_noise*log(1-est_Ki_Two_Comp_noise/est_F_noise);

% Create the estimated filter out of the non-linear parameters estimation
Filter_estimation_result  = est_F_noise * Larsson_Filter( time_vec_minutes', est_F_noise , ...
    est_Vb_Two_Comp_noise, E_Two_Comp_est , Ve_Two_Comp_est, Hct(iter_num));

if (plot_flag)
    % Check to see if one of the filters explodes and will ruin visualization
    All_filters           = zeros(7,max(size(larss_filter)));
    All_filters(1,:)      = larss_filter;
    All_filters(2,:)      = est_larss_filter_Wiener_noise;
    All_filters(3,:)      = ridge_regression_larss_result;
    All_filters(4,:)      = b_spline_larss_result;
    All_filters(5,:)      = b_spline_larss_result_1st_deriv;
    All_filters(6,:)      = b_spline_larss_result_2nd_deriv;
    All_filters(7,:)      = b_PCA_larss_result_2nd_deriv;
    Max_values            = max(All_filters,[],2);
    Max_original          = max(larss_filter);
    Bad_indices           = find(Max_values > 5*Max_original);
    % Zero exloding filters so they won't ruin the plot
    All_filters(Bad_indices, :) = 0;
    
    fig_num = figure;
    subplot(2,1,1);
    hold on;
    h1  = plot(time_vec_minutes,All_filters(1,:),'g');
    h2  = plot(time_vec_minutes,All_filters(1,:),'g*');
    h3  = plot(time_vec_minutes,All_filters(2,:),'b');
    h4  = plot(time_vec_minutes,All_filters(2,:),'b+');
    h5  = plot(time_vec_minutes,All_filters(3,:),'k');
    h6  = plot(time_vec_minutes,All_filters(3,:),'ko');
    h7  = plot(time_vec_minutes,All_filters(4,:),'c');
    h8  = plot(time_vec_minutes,All_filters(4,:),'cx');
    h9  = plot(time_vec_minutes,All_filters(5,:),'m');
    h10 = plot(time_vec_minutes,All_filters(5,:),'md');
    h11 = plot(time_vec_minutes,All_filters(6,:) ,'r');
    h12 = plot(time_vec_minutes,All_filters(6,:) ,'rs');
    h13 = plot(time_vec_minutes,All_filters(7,:) ,'y');
    h14 = plot(time_vec_minutes,All_filters(7,:) ,'yh');
    
    % Warn in the title in case one filter exploded
    if (isempty(Bad_indices))
        title('Original and estimated h(t) - Larsson','FontWeight','bold');
    else
        title('Original and estimated h(t) - Larsson - REMOVED BAD FILTERS','FontWeight','bold');
    end
    
    xlabel('Time [Min]');
    hold off;
    legend([h2 h4 h6 h8 h10 h12 h14],'Orig. h(t)',...
        'Est. h(t) - Wiener','Est. h(t) - Ridge','Est. h(t) - Spline','Est. h(t) - Spline 1st deriv','Est. h(t) - Spline 2nd deriv','Est. h(t) - PCA 2nd deriv');
    
    subplot(2,1,2);
    hold on;
    
    AIF_filtered_by_est_larrson  = filter(est_larss_filter_Wiener_noise*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num));
    AIF_filtered_by_est_ridge    = filter(ridge_regression_larss_result*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num));
    AIF_filtered_by_est_spline_1 = filter(b_spline_larss_result*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num));
    AIF_filtered_by_est_spline_2 = filter(b_spline_larss_result_1st_deriv*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num));
    AIF_filtered_by_est_spline_3 = filter(b_spline_larss_result_2nd_deriv*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num));
    AIF_filtered_by_est_spline_4 = filter(b_PCA_larss_result_2nd_deriv*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num));
    
    h1  = plot(time_vec_minutes,Sim_Ct_larss_kernel_noise,'g');
    h2  = plot(time_vec_minutes,Sim_Ct_larss_kernel_noise,'g*');
    h3  = plot(time_vec_minutes,AIF_filtered_by_est_larrson,'b');
    h4  = plot(time_vec_minutes,AIF_filtered_by_est_larrson,'b+');
    h5  = plot(time_vec_minutes,AIF_filtered_by_est_ridge,'k');
    h6  = plot(time_vec_minutes,AIF_filtered_by_est_ridge,'ko');
    h7  = plot(time_vec_minutes,AIF_filtered_by_est_spline_1,'c');
    h8  = plot(time_vec_minutes,AIF_filtered_by_est_spline_1,'cx');
    h9  = plot(time_vec_minutes,AIF_filtered_by_est_spline_2,'m');
    h10 = plot(time_vec_minutes,AIF_filtered_by_est_spline_2,'md');
    h11 = plot(time_vec_minutes,AIF_filtered_by_est_spline_3,'r');
    h12 = plot(time_vec_minutes,AIF_filtered_by_est_spline_3,'rs');
    h13 = plot(time_vec_minutes,AIF_filtered_by_est_spline_4,'y');
    h14 = plot(time_vec_minutes,AIF_filtered_by_est_spline_4,'yh');
    
    title('Original and estimated Ct(t)','FontWeight','bold');
    xlabel('Time [Min]');
    hold off;
    legend([h2 h4 h6 h8 h10 h12 h14],'Orig. h(t)',...
        'Est. Ct(t) - Wiener','Est. Ct(t) - Ridge','Est. Ct(t) - Spline','Est. Ct(t) - Spline 1st deriv','Est. Ct(t) - Spline 2nd deriv','Est. Ct(t) - PCA 2nd deriv');
    
    % Print result to PDF
    [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Filters_Larss.png', './Run_Output/',...
        'Estimated Filters after de-convolution - Larsson', 'EstFiltersLarss');
    
    % Plot the estimated filter using the non-linear paramters estimation
    fig_num = figure;
    hold on;
    h1 = plot(time_vec_minutes,larss_filter,'g');
    h2 = plot(time_vec_minutes,larss_filter,'g*');
    h3 = plot(time_vec_minutes,Final_Filter_Estimation_Larss ,'r');
    h4 = plot(time_vec_minutes,Final_Filter_Estimation_Larss ,'rs');
    h5 = plot(time_vec_minutes,Filter_estimation_result,'m');
    h6 = plot(time_vec_minutes,Filter_estimation_result,'md');
    hold off;
    
    title('Estimated h(t) and non-linear fit','FontWeight','bold');
    xlabel('Time [Min]');
    hold off;
    legend([h2 h4 h6],'Original h(t)','Est. h(t)','Fitted h(t) - Non-Linear');
    
    % Print result to PDF
    [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Filter_and_Fitted_one.png', './Run_Output/',...
        'Estimated Filter and non-linear fit', 'EstFiltersLarss');
    
    % Display Parameters estimation result vs original
    fig_num = figure;
    
    annotation('textbox',...
        [0 0 1 1],...
        'String',{'Larsson Filter Original Params:','',['F         = ' num2str(F,'%.2f') '        [mL/100g/min]'],...
        ['Delay  = ' num2str(additional_AIF_delay_sec,'%.2f') '          [sec]'],...
        ['Ki        = ' num2str(Ki,'%.2f') '          [mL/100g/min]'],...
        ['PS      = ' num2str(PS,'%.2f') '           [mL/100g/min]'],...
        ['Vb       = ' num2str(Vb_larss,'%.2f') '        [mL/100g]'],...
        ['Ve       = ' num2str(Ve_larss,'%.2f') '          [mL/100g]'],...
        ['Vd       = ' num2str(Vd,'%.2f') '        [mL/100g]'],...
        ['E         = ' num2str(E,'%.2f')],...
        ['MTT   = ' num2str(MTT,'%.2f') '           [sec]'],...
        ['SNR   = ' num2str(SNR_single,'%.2f')],'',...
        ...
        'Estimated Params:','',...
        ['F est                       = ' num2str(est_F_noise,'%.2f') '      [mL/100g/min]'],...
        ['Delay                      = ' num2str(est_Delay_sec_noise,'%.2f') ...
        '       [Sec]'],...
        ['Delay-Gauss          = ' num2str(est_t_d_noise,'%.2f') ...
        '       [Sec]'],...
        ['Ki-Patlak                 = ' num2str(est_Ki_Patlak_noise,'%.2f') '       [mL/100g/min]'],...
        ['Ki-Two-Comp         = ' num2str(est_Ki_Two_Comp_noise,'%.2f') '        [mL/100g/min]'],...
        ['E-Patlak                  = ' num2str(E_Patlak_est,'%.2f')],...
        ['E-Two-Comp          = ' num2str(E_Two_Comp_est,'%.2f')],...
        ['PS-est                     = ' num2str(est_PS_noise,'%.2f') '       [mL/100g/min]'],...
        ['Vb-Patlak                = ' num2str( est_Vb_Patlak_noise,'%.2f') '     [mL/100g]'],...
        ['Vb-Two-Comp         = ' num2str(est_Vb_Two_Comp_noise,'%.2f') '     [mL/100g]'],...
        ['Vd-est                     = ' num2str(est_Vd_noise,'%.2f') '      [mL/100g]'],...
        ['Vd-normal-tis-est     = ' num2str(est_Vd_normal_tis_noise,'%.2f') '       [mL/100g]'],...
        ['MTT-est                  = ' num2str(est_MTT_noise,'%.2f') '          [Sec]'],...
        ['MTT-normal-tis-est = ' num2str(est_MTT_normal_tis_noise,'%.2f') '          [Sec]']},...
        'FontSize',8,...
        'FontName','Arial',...
        'LineStyle','-',...
        'EdgeColor',[0 0 0],...
        'LineWidth',2,...
        'BackgroundColor',[1 1 1],...
        'Color',[0.0 0.0 0]);
    
    % Print result to PDF
    [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Estimaed_vs_Real_Larsson_Params.png', './Run_Output/',...
        'Estimated Vs. Real Params - Larsson', 'ParamsResults');
    
end

% Return vals
Return_Struct.est_F_noise                        = est_F_noise;
Return_Struct.est_Vb_Patlak_noise                = est_Vb_Patlak_noise;
Return_Struct.E_Patlak_est                       = E_Patlak_est;
Return_Struct.est_Ki_Patlak_noise                = est_Ki_Patlak_noise;
Return_Struct.est_Delay_sec_using_Gaussian_noise = est_Delay_sec_using_Gaussian_noise;
Return_Struct.est_Delay_sec_noise                = est_Delay_sec_noise;
Return_Struct.est_Ki_Two_Comp_noise              = est_Ki_Two_Comp_noise;
Return_Struct.est_E_noise                        = est_E_noise;
Return_Struct.est_PS_noise                       = est_PS_noise;
Return_Struct.est_Vb_Two_Comp_noise              = est_Vb_Two_Comp_noise;
Return_Struct.est_Vd_noise                       = est_Vd_noise;
Return_Struct.est_Vd_normal_tis_noise            = est_Vd_normal_tis_noise;
Return_Struct.est_MTT_noise                      = est_MTT_noise;
Return_Struct.est_MTT_normal_tis_noise           = est_MTT_normal_tis_noise;

if Check_Sourbron_Estimate
    Return_Struct.est_F_Two_Comp_Sourbron_noise  = est_F_Two_Comp_Sourbron_noise;
    Return_Struct.est_Vb_Two_Comp_Sourbron_noise = est_Vb_Two_Comp_Sourbron_noise; 
    Return_Struct.est_Ki_Two_Comp_Sourbron_noise = est_Ki_Two_Comp_Sourbron_noise;
    Return_Struct.Ve_Two_Comp_Sourbron_est       = Ve_Two_Comp_Sourbron_est;   
end

end