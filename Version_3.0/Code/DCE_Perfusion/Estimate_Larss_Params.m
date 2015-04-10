function [ Return_Struct, idx_fig ] = Estimate_Larss_Params( Sim_Struct, ht_Struct, est_t_d_noise, Verbosity, iter_num, avg_num, idx_fig)

% Take from struct variables used in local function
Filter_Est_Chosen               = Sim_Struct.Filter_Est_Chosen;
FMS_Algorithm                   = Sim_Struct.FMS_Algorithm;
Init_Guess_Gaussian             = Sim_Struct.Init_Guess_Gaussian;
LowerBound_Gauss                = Sim_Struct.LowerBound_Gauss;
UpperBound_Gauss                = Sim_Struct.UpperBound_Gauss;
Sim_AIF_delayed_no_noise        = Sim_Struct.Sim_AIF_delayed_no_noise;
Ktrans                          = Sim_Struct.Ktrans;
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
Random_init_F_guess_4_Sourbron  = Sim_Struct.Random_init_F_guess_4_Sourbron;
adjusted_larsson                = Sim_Struct.Adjusted_Larsson_Model;
F_low                           = Sim_Struct.F_low;
F_max                           = Sim_Struct.F_max;
Vb_low                          = Sim_Struct.Vb_low;
Vb_max                          = Sim_Struct.Vb_max;
Ve_low                          = Sim_Struct.Ve_low;
Ve_max                          = Sim_Struct.Ve_max;
E_low                           = Sim_Struct.E_low;
E_max                           = Sim_Struct.E_max;
init_Ve_guess                   = Sim_Struct.init_Ve_guess;

est_larss_filter_Wiener_noise   = Sim_Struct.est_larss_filter_Wiener_noise;
ridge_regression_larss_result   = ht_Struct.ridge_regression_larss_result;
b_spline_larss_result           = ht_Struct.b_spline_larss_result;
b_spline_larss_result_1st_deriv = ht_Struct.b_spline_larss_result_1st_deriv;
b_spline_larss_result_2nd_deriv = ht_Struct.b_spline_larss_result_2nd_deriv;
b_PCA_larss_result              = ht_Struct.b_PCA_larss_result;
b_PCA_larss_result_1st_deriv    = ht_Struct.b_PCA_larss_result_1st_deriv;
b_PCA_larss_result_2nd_deriv    = ht_Struct.b_PCA_larss_result_2nd_deriv;

Final_Filter_Estimation_Larss   = ht_Struct.Final_Filter_Estimation_Larss;

if strcmp(Verbosity,'Full')
    display('-I- Estimating Larsson parameters...');
end


% F estimated is the maximum value of F*IRF
est_F_noise         = max(Final_Filter_Estimation_Larss);

% Delay of the AIF will be calculated according to the place of the
% maximum value of F*IRF
max_index           = Final_Filter_Estimation_Larss == est_F_noise;

% Translate from minutes
if (Correct_estimation_due_to_delay)
    est_Delay_sec_noise =  est_delay_by_AIF_correct;
else
    est_Delay_sec_noise =  time_vec_minutes(max_index) * 60;
end


%% Patalk Estimation
[est_Ktrans_Patlak_noise, est_Vb_Patlak_noise ,est_E_Patlak_noise, est_MTT_Patlak_noise, idx_fig] = Patlak_Estimation(Sim_Struct, Sim_AIF_with_noise, Sim_Ct_larss_kernel_noise, est_F_noise, Verbosity, iter_num, avg_num, idx_fig);

%% Using Murase estimation for Vb,Ktrans - Alternative to Patlak (Using estimated F)

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

% -------------- Ktrans estimation ---------------
Y_vec_Ktrans            = ( Vb_murase_estimate*est_F_noise*cumtrapz(time_vec_minutes,Sim_AIF_delayed_no_noise) - ...
    - Sim_Ct_larss_kernel - est_F_noise* cumtrapz(time_vec_minutes,Sim_Ct_larss_kernel) ) ./ ...
    ( cumtrapz(time_vec_minutes,Sim_Ct_larss_kernel));

Ktrans_murase_estimate  = mean(Y_vec_Ktrans(stable_idx));


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
    h1 = plot(time_vec_minutes,Y_vec_Ktrans,'k*');
    h2 = plot(time_vec_minutes(stable_idx),Y_vec_Ktrans(stable_idx),'og');
    hold off;
    title('Ktrans Murase Estimation','FontWeight','bold');
    xlabel('Time [Min]');
    legend([h1 h2],'Ktrans vector','Stable points for estimation');
    
    % Display Larsson's parameters
    annotation('textbox',...
        [0.7 0.09 0.3 0.1],...
        'String',{['Ktrans = ' num2str(Ktrans) '  [mL/100g]'],...
        ['Ktrans Murase = ' num2str(Ktrans_murase_estimate,'%.2f') '  [mL/100g]'],...
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
    [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Vb_Ktrans_Murase_Estimate.png', './Run_Output/',...
        'Murase Estimate for Vb, Ktrans', 'VbKtransMuraseEst');
    
end

%% Estimate using two compartment model using initial Patlak's guess

% Initial Guess for non-linear curve fitting for Larsson (Vb, E, Ve)
Init_Guess_Larsson = [ est_Vb_Patlak_noise est_E_Patlak_noise init_Ve_guess];
%Init_Guess_Larsson = [50 0.15 10];

% The analytic funcational of a Larsson function
if (adjusted_larsson)
    Larsson_function             = @(x,t) Adjusted_Larsson_Filter( t, est_F_noise, x(1), x(2), x(3));
    Larsson_function_4_Sourbron  = @(x,t) Adjusted_Larsson_Filter( t, x(1)       , x(2), x(3), x(4));
else
    Larsson_function             = @(x,t) Larsson_Filter( t, est_F_noise, x(1), x(2), x(3), Hct(iter_num));
    Larsson_function_4_Sourbron  = @(x,t) Larsson_Filter( t, x(1)       , x(2), x(3), x(4), Hct(iter_num));
end

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
        if (Random_init_F_guess_4_Sourbron)
            F_init_guess        = rand*100;
        else
            F_init_guess        = est_F_noise;
        end
        
        Init_Guess_Sourbron = [F_init_guess Init_Guess_Larsson];
        
        LowerBound_Sourbron = [F_low LowerBound_Larsson];
        UpperBound_Sourbron = [F_max UpperBound_Larsson];
        
        % F, Vb, E, Ve
        [est_params_Sourbron_noise,~,~,~,~] = ...
            lsqcurvefit(Larsson_function_4_Sourbron,Init_Guess_Sourbron,time_vec_minutes',Final_Filter_Estimation_Larss/F_init_guess,...
            LowerBound_Sourbron,UpperBound_Sourbron,algorithm_options);
        
        
    end
    
elseif strcmp(FMS_Algorithm,'levenberg-marquardt')
    
    %     Bounded_Larsson_function             = @(x) BoundFunc(Larsson_function,x,LowerBound_Larsson,UpperBound_Larsson, time_vec_minutes, est_F_noise);
    %     [~, Init_Guess_Larsson_Transformed]  = BoundFunc(Larsson_function,Init_Guess_Larsson,LowerBound_Larsson,UpperBound_Larsson, time_vec_minutes, est_F_noise);
    %
    %     Larsson_function             = @(x,t) Adjusted_Larsson_Filter( t, est_F_noise, x(1), x(2), x(3));
    %
    %     ydata    = (Final_Filter_Estimation_Larss'/est_F_noise);
    %     RMSCost  = @(vec) sum(vec.^2);
    %     CostFunc = @(x) RMSCost( Bounded_Larsson_function(x) - ydata );
    %     best_transformedParams    = fminsearch(CostFunc, Init_Guess_Larsson_Transformed);
    %     % Transform back
    %     %est_params_Larsson_noise = ReverseTransformation(best_transformedParams, LowerBound_Larsson, UpperBound_Larsson);
    
    Unbounded_Larsson_function = @(x,t) Larsson_function(BoundFunc(x,LowerBound_Larsson,UpperBound_Larsson),t);
    
    [est_params_Larsson_noise,residue_norm_Larsson_noise,residual_Larsson_noise,exitflag_Larsson_noise,algo_info_Larsson_noise] = ...
        lsqcurvefit(Unbounded_Larsson_function,Init_Guess_Larsson,time_vec_minutes',Final_Filter_Estimation_Larss/est_F_noise,...
        [],[],algorithm_options);
    
    est_params_Larsson_noise = BoundFunc(est_params_Larsson_noise,LowerBound_Larsson,UpperBound_Larsson);
    
elseif strcmp(FMS_Algorithm,'fminsearch')
    
    Unbounded_Larsson_function = @(x,t) Larsson_function(BoundFunc(x,LowerBound_Larsson,UpperBound_Larsson),t);
    ydata    = Final_Filter_Estimation_Larss/est_F_noise;
    RMSCost  = @(vec) sum(vec.^2);
    CostFunc = @(x) RMSCost( Unbounded_Larsson_function(x,time_vec_minutes') - ydata );
    best_transformedParams    = fminsearch(CostFunc, Init_Guess_Larsson, Sim_Struct.algorithm_options);
    est_params_Larsson_noise = BoundFunc(best_transformedParams,LowerBound_Larsson,UpperBound_Larsson);
    
    
    
    %     [est_params_Larsson_noise,residue_norm_Larsson_noise,residual_Larsson_noise,exitflag_Larsson_noise,algo_info_Larsson_noise] = ...
    %         lsqcurvefit(Larsson_function,Init_Guess_Larsson,time_vec_minutes',Final_Filter_Estimation_Larss/est_F_noise,...
    %         [],[],algorithm_options);
    
    % Estimate gaussian parameters when Larss filter is used
    if ~Ignore_Gaussian_Calculation
        
        [est_params_Larss_Usiung_Gauss_noise,~,~,~,~] = ...
            lsqcurvefit(Gaussian_function,Init_Guess_Gaussian,time_vec_minutes',Final_Filter_Estimation_Larss,[],[],algorithm_options);
        
    end
    
    if Check_Sourbron_Estimate
        
        % Gues initial F randomly
        if (Random_init_F_guess_4_Sourbron)
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
est_E_Two_Comp_noise            = est_params_Larsson_noise(2);
est_Ktrans_Two_Comp_noise       = est_E_Two_Comp_noise * est_F_noise;
est_Ve_Two_Comp_noise           = est_params_Larsson_noise(3);

% Assigning two compartment parameters estimation using Sourbron
if Check_Sourbron_Estimate
    est_F_Two_Comp_Sourbron_noise               = est_params_Sourbron_noise(1);
    est_Vb_Two_Comp_Sourbron_noise              = est_params_Sourbron_noise(2);
    E_Two_Comp_Sourbron_est                     = est_params_Sourbron_noise(3);
    est_Ktrans_Two_Comp_Sourbron_noise              = E_Two_Comp_Sourbron_est * est_F_Two_Comp_Sourbron_noise;
    Ve_Two_Comp_Sourbron_est                    = est_params_Sourbron_noise(4);
end

% Estimate MTT (also in normal tissue)
est_IRF                         = Final_Filter_Estimation_Larss / est_F_noise;
est_MTT_noise                   = cumtrapz(time_vec_minutes,est_IRF);
est_MTT_noise                   = est_MTT_noise(end);
est_MTT_normal_tis_noise        = est_Vb_Patlak_noise / est_F_noise;
% Estimate Vd (also in normal tissue)
est_Vd_noise                    = est_F_noise * est_MTT_noise;
est_Vd_normal_tis_noise         = est_F_noise * est_MTT_normal_tis_noise;
% Estimate E
est_E_noise                     = est_E_Two_Comp_noise;
% Estimate PS
est_PS_noise                    = -est_F_noise*log(1-est_Ktrans_Two_Comp_noise/est_F_noise);

% Create the estimated filter out of the non-linear parameters estimation
if (adjusted_larsson)
    Filter_estimation_result  = est_F_noise * Adjusted_Larsson_Filter( time_vec_minutes', est_F_noise , ...
        est_Vb_Two_Comp_noise, est_E_Two_Comp_noise , est_Ve_Two_Comp_noise);
else
    Filter_estimation_result  = est_F_noise * Larsson_Filter( time_vec_minutes', est_F_noise , ...
        est_Vb_Two_Comp_noise, est_E_Two_Comp_noise , est_Ve_Two_Comp_noise, Hct(iter_num));
end

if (plot_flag)
    % Check to see if one of the filters explodes and will ruin visualization
    All_filters           = zeros(9,max(size(larss_filter)));
    All_filters(1,:)      = larss_filter;
    All_filters(2,:)      = est_larss_filter_Wiener_noise;
    All_filters(3,:)      = ridge_regression_larss_result;
    All_filters(4,:)      = b_spline_larss_result;
    All_filters(5,:)      = b_spline_larss_result_1st_deriv;
    All_filters(6,:)      = b_spline_larss_result_2nd_deriv;
    All_filters(7,:)      = b_PCA_larss_result;
    All_filters(8,:)      = b_PCA_larss_result_1st_deriv;
    All_filters(9,:)      = b_PCA_larss_result_2nd_deriv;
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
    h13 = plot(time_vec_minutes,All_filters(7,:) ,'b--');
    h14 = plot(time_vec_minutes,All_filters(7,:) ,'bh');
    h15 = plot(time_vec_minutes,All_filters(8,:) ,'k--');
    h16 = plot(time_vec_minutes,All_filters(8,:) ,'kh');
    h17 = plot(time_vec_minutes,All_filters(9,:) ,'y');
    h18 = plot(time_vec_minutes,All_filters(9,:) ,'yh');
    
    % Warn in the title in case one filter exploded
    if (isempty(Bad_indices))
        title('True and estimated h(t) - Larsson','FontWeight','bold');
    else
        title('True and estimated h(t) - Larsson - REMOVED BAD FILTERS','FontWeight','bold');
    end
    
    xlabel('Time [Min]');
    hold off;
    legend([h2 h4 h6 h8 h10 h12 h14 h16 h18],'Orig. h(t)',...
        'Est. h(t) - Wiener','Est. h(t) - Ridge','Est. h(t) - Spline','Est. h(t) - Spline 1st deriv','Est. h(t) - Spline 2nd deriv','Est. h(t) - PCA','Est. h(t) - PCA 1st deriv','Est. h(t) - PCA 2nd deriv');
    
    subplot(2,1,2);
    hold on;
    if Sim_Struct.ignore_time_delta
        AIF_filtered_by_est_larrson          = filter(est_larss_filter_Wiener_noise,1,Sim_AIF_with_noise(:,iter_num,avg_num));
        AIF_filtered_by_est_ridge            = filter(ridge_regression_larss_result,1,Sim_AIF_with_noise(:,iter_num,avg_num));
        AIF_filtered_by_est_spline_no_deriv  = filter(b_spline_larss_result,1,Sim_AIF_with_noise(:,iter_num,avg_num));
        AIF_filtered_by_est_spline_1st_deriv = filter(b_spline_larss_result_1st_deriv,1,Sim_AIF_with_noise(:,iter_num,avg_num));
        AIF_filtered_by_est_spline_2nd_deriv = filter(b_spline_larss_result_2nd_deriv,1,Sim_AIF_with_noise(:,iter_num,avg_num));
        AIF_filtered_by_est_PCA              = filter(b_PCA_larss_result,1,Sim_AIF_with_noise(:,iter_num,avg_num));
        AIF_filtered_by_est_PCA_1st_deriv    = filter(b_PCA_larss_result_1st_deriv,1,Sim_AIF_with_noise(:,iter_num,avg_num));
        AIF_filtered_by_est_PCA_2nd_deriv    = filter(b_PCA_larss_result_2nd_deriv,1,Sim_AIF_with_noise(:,iter_num,avg_num));
    else
        AIF_filtered_by_est_larrson          = filter(est_larss_filter_Wiener_noise*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num));
        AIF_filtered_by_est_ridge            = filter(ridge_regression_larss_result*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num));
        AIF_filtered_by_est_spline_no_deriv  = filter(b_spline_larss_result*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num));
        AIF_filtered_by_est_spline_1st_deriv = filter(b_spline_larss_result_1st_deriv*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num));
        AIF_filtered_by_est_spline_2nd_deriv = filter(b_spline_larss_result_2nd_deriv*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num));
        AIF_filtered_by_est_PCA              = filter(b_PCA_larss_result*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num));
        AIF_filtered_by_est_PCA_1st_deriv    = filter(b_PCA_larss_result_1st_deriv*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num));
        AIF_filtered_by_est_PCA_2nd_deriv    = filter(b_PCA_larss_result_2nd_deriv*min_interval,1,Sim_AIF_with_noise(:,iter_num,avg_num));
    end
    
    h1  = plot(time_vec_minutes,Sim_Ct_larss_kernel_noise,'g');
    h2  = plot(time_vec_minutes,Sim_Ct_larss_kernel_noise,'g*');
    h3  = plot(time_vec_minutes,AIF_filtered_by_est_larrson,'b');
    h4  = plot(time_vec_minutes,AIF_filtered_by_est_larrson,'b+');
    h5  = plot(time_vec_minutes,AIF_filtered_by_est_ridge,'k');
    h6  = plot(time_vec_minutes,AIF_filtered_by_est_ridge,'ko');
    h7  = plot(time_vec_minutes,AIF_filtered_by_est_spline_no_deriv,'c');
    h8  = plot(time_vec_minutes,AIF_filtered_by_est_spline_no_deriv,'cx');
    h9  = plot(time_vec_minutes,AIF_filtered_by_est_spline_1st_deriv,'m');
    h10 = plot(time_vec_minutes,AIF_filtered_by_est_spline_1st_deriv,'md');
    h11 = plot(time_vec_minutes,AIF_filtered_by_est_spline_2nd_deriv,'r');
    h12 = plot(time_vec_minutes,AIF_filtered_by_est_spline_2nd_deriv,'rs');
    h13 = plot(time_vec_minutes,AIF_filtered_by_est_PCA ,'b--');
    h14 = plot(time_vec_minutes,AIF_filtered_by_est_PCA ,'bh');
    h15 = plot(time_vec_minutes,AIF_filtered_by_est_PCA_1st_deriv,'k--');
    h16 = plot(time_vec_minutes,AIF_filtered_by_est_PCA_1st_deriv,'kh');
    h17 = plot(time_vec_minutes,AIF_filtered_by_est_PCA_2nd_deriv,'y');
    h18 = plot(time_vec_minutes,AIF_filtered_by_est_PCA_2nd_deriv,'yh');
    
    title('True and estimated Ct(t)','FontWeight','bold');
    xlabel('Time [Min]');
    hold off;
    legend([h2 h4 h6 h8 h10 h12 h14 h16 h18],'Orig. h(t)',...
        'Est. Ct(t) - Wiener','Est. Ct(t) - Ridge','Est. Ct(t) - Spline','Est. Ct(t) - Spline 1st deriv','Est. Ct(t) - Spline 2nd deriv','Est. Ct(t) - PCA','Est. Ct(t) - PCA 1st deriv','Est. Ct(t) - PCA 2nd deriv');
    
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
    legend([h2 h4 h6],'True h(t)','Est. h(t)','Fitted h(t) - Non-Linear');
    
    % Print result to PDF
    [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Filter_and_Fitted_one.png', './Run_Output/',...
        'Estimated Filter and non-linear fit', 'EstFiltersLarss');
    
    % Display Parameters estimation result vs original
    fig_num = figure;
    
    annotation('textbox',...
        [0 0 1 1],...
        'String',{'Larsson Filter True Params:','',['F         = ' num2str(F,'%.2f') '        [mL/100g/min]'],...
        ['Delay  = ' num2str(additional_AIF_delay_sec,'%.2f') '          [sec]'],...
        ['Ktrans        = ' num2str(Ktrans,'%.2f') '          [mL/100g/min]'],...
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
        ['Ktrans-Patlak                 = ' num2str(est_Ktrans_Patlak_noise,'%.2f') '       [mL/100g/min]'],...
        ['Ktrans-Two-Comp         = ' num2str(est_Ktrans_Two_Comp_noise,'%.2f') '        [mL/100g/min]'],...
        ['E-Patlak                  = ' num2str(est_E_Patlak_noise,'%.2f')],...
        ['E-Two-Comp          = ' num2str(est_E_Two_Comp_noise,'%.2f')],...
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
Return_Struct.E_Patlak_est                       = est_E_Patlak_noise;
Return_Struct.est_Ktrans_Patlak_noise                = est_Ktrans_Patlak_noise;
Return_Struct.est_Delay_sec_using_Gaussian_noise = est_Delay_sec_using_Gaussian_noise;
Return_Struct.est_Delay_sec_noise                = est_Delay_sec_noise;
Return_Struct.est_Ktrans_Two_Comp_noise              = est_Ktrans_Two_Comp_noise;
Return_Struct.est_E_noise                        = est_E_noise;
Return_Struct.est_PS_noise                       = est_PS_noise;
Return_Struct.est_Vb_Two_Comp_noise              = est_Vb_Two_Comp_noise;
Return_Struct.est_Ve_Two_Comp_noise              = est_Ve_Two_Comp_noise;
Return_Struct.est_Vd_noise                       = est_Vd_noise;
Return_Struct.est_Vd_normal_tis_noise            = est_Vd_normal_tis_noise;
Return_Struct.est_MTT_noise                      = est_MTT_noise;
Return_Struct.est_MTT_normal_tis_noise           = est_MTT_normal_tis_noise;

if Check_Sourbron_Estimate
    Return_Struct.est_F_Two_Comp_Sourbron_noise  = est_F_Two_Comp_Sourbron_noise;
    Return_Struct.est_Vb_Two_Comp_Sourbron_noise = est_Vb_Two_Comp_Sourbron_noise;
    Return_Struct.est_Ktrans_Two_Comp_Sourbron_noise = est_Ktrans_Two_Comp_Sourbron_noise;
    Return_Struct.Ve_Two_Comp_Sourbron_est       = Ve_Two_Comp_Sourbron_est;
end

end