function [ est_delay_by_AIF_correct, Sim_AIF_with_noise_Regul_shifted, Final_Filter_Estimation_Larss, idx_fig ] = AIF_Delay_Correct( Sim_Struct, ht_Struct, Final_Filter_Estimation_Larss, Sim_Ct_larss_Regul_noise,Verbosity, iter_num, avg_num, idx_fig)

% Take from struct variables used in local function
normalize                       = Sim_Struct.normalize;
B_mat                           = Sim_Struct.B_mat;
B_PCA                           = Sim_Struct.PCA_B_mat;
plot_L_Curve                    = Sim_Struct.plot_L_Curve;
Derivative_Time_Devision        = Sim_Struct.Derivative_Time_Devision;
lambda_vec_larss                = Sim_Struct.lambda_vec_larss;
min_interval                    = Sim_Struct.min_interval(iter_num);
time_vec_minutes                = double(Sim_Struct.time_vec_minutes);
Upsampling_resolution           = Sim_Struct.Upsampling_resolution;
Max_Time_Delay                  = Sim_Struct.Max_Time_Delay;
Min_Time_Delay                  = Sim_Struct.Min_Time_Delay;
Use_Upsampling_Delay_Comp       = Sim_Struct.Use_Upsampling_Delay_Comp;
LowerBound_Larsson              = Sim_Struct.LowerBound_Larsson;
UpperBound_Larsson              = Sim_Struct.UpperBound_Larsson;
algorithm_options               = Sim_Struct.algorithm_options;
Hct                             = Sim_Struct.Hct;
RMS_Smooth_Around_Bolus         = Sim_Struct.RMS_Smooth_Around_Bolus;
RMS_Smooth                      = Sim_Struct.RMS_Smooth;
Diff_From_Bolus                 = Sim_Struct.Diff_From_Bolus;
additional_AIF_delay_sec        = Sim_Struct.additional_AIF_delay_sec(iter_num);
BiExp2CTC_RMS_Ratio             = Sim_Struct.BiExp2CTC_RMS_Ratio;
plot_flag                       = Sim_Struct.plot_flag;
adjusted_larsson                = Sim_Struct.Adjusted_Larsson_Model;
Filter_Est_Chosen               = Sim_Struct.Filter_Est_Chosen;
RealData_Flag                   = Sim_Struct.RealData_Flag;
Simple_AIF_Delay_Correct        = Sim_Struct.Simple_AIF_Delay_Correct;
init_Ve_guess                   = Sim_Struct.init_Ve_guess;
FMS_Algorithm                   = Sim_Struct.FMS_Algorithm;

%Final_Filter_Estimation_Larss   = ht_Struct.Final_Filter_Estimation_Larss;
%Sim_Ct_larss_Regul_noise        = ht_Struct.Sim_Ct_larss_Regul_noise;
Sim_AIF_with_noise_Regul        = ht_Struct.Sim_AIF_with_noise_Regul;
Sim_Ct_larss_Regul              = ht_Struct.Sim_Ct_larss_Regul;
Conv_X_no_noise                 = ht_Struct.Conv_X_no_noise;

% ---------------------------------------------------------------------
%                 Regular correction by looking on h(t) shift
% ---------------------------------------------------------------------
est_simple_delay       = NaN; % Initiate with NaN

if Simple_AIF_Delay_Correct
    
    [~, max_idx] = max(Final_Filter_Estimation_Larss);
    
    % If the maximum is not the first value, we might think there is a delay
    if (max_idx ~= 1)
        
        [peak_val, peak_idx]       = findpeaks(Final_Filter_Estimation_Larss);
        
        % We expect a peak in case of delay
        if ~isempty(peak_idx)
            
            % The estimated delay in minutes
            est_simple_delay                   = time_vec_minutes(peak_idx(1));
            
            % Shift the AIF according to estimation
            shift_index                                             = peak_idx(1);
            Sim_AIF_with_noise_Regul_shifted                        = zeros(size(Sim_AIF_with_noise_Regul));
            Sim_AIF_with_noise_Regul_shifted(1:end-shift_index)     = Sim_AIF_with_noise_Regul(shift_index+1:end);
            Sim_AIF_with_noise_Regul_shifted(end-shift_index+1:end) = Sim_AIF_with_noise_Regul_shifted(end-shift_index); % Duplicate the last values?
            
            % Create new convolution matrix for
            [ Conv_X_shifted ] = Convolution_Matrix( min_interval, Sim_AIF_with_noise_Regul_shifted );
            
            if RealData_Flag
                [ ~, ~, ~, Final_Filter_Estimation_Larss, ~, ~, ~, idx_fig ] =  ...
                    Regularization_Methods_Simulation( Sim_Ct_larss_Regul, Sim_Ct_larss_Regul_noise, Conv_Matrix, Conv_Matrix_no_noise, time_vec_minutes, lambda_vec_larss, normalize, min_interval, B_mat, B_PCA, plot_L_Curve, idx_fig, filter_type, Derivative_Time_Devision, plot_flag, RealData_Flag );
            else
                
                [ridge_regression_result, b_spline_result, b_spline_result_1st_deriv, b_spline_result_2nd_deriv, b_PCA_result_1st_deriv, b_PCA_result_2nd_deriv, idx_fig]...
                    = Regularization_Methods_Simulation(Sim_Ct_larss_Regul, Sim_Ct_larss_Regul_noise,Conv_X_shifted,Conv_X_no_noise,time_vec_minutes,...
                    lambda_vec_larss, normalize, min_interval, B_mat, plot_L_Curve, idx_fig , 'Larss' , Derivative_Time_Devision, plot_flag, RealData_Flag );
                
                [ Final_Filter_Estimation_Larss ] = Choose_Needed_Ht( Filter_Est_Chosen, est_larss_filter_Wiener_noise, ridge_regression_result, b_spline_result, b_spline_result_1st_deriv, b_spline_result_2nd_deriv, b_PCA_larss_result, b_PCA_result_1st_deriv, b_PCA_result_2nd_deriv);
            end
            
            
        end
    else
        est_simple_delay = 0; % If the maximum is the first delay, there is no delay
    end
    
    est_delay_by_AIF_correct = est_simple_delay;
    
    return;
    
    % ---------------------------------------------------------------------
    %                 Correction by setting multiple optional time delays
    % ---------------------------------------------------------------------
else
    
    %     [peak_val, peak_idx]       = findpeaks(b_spline_larss_result_2nd_deriv);
    %     est_delay_by_spline_result = time_vec_minutes(peak_idx(1));
    %     figure;
    %     plot(time_vec_minutes*60,b_spline_larss_result_2nd_deriv);
    %     xlabel('Time [Sec]');
    %     title('Estimated Filter');
    %
    %     % Check min RMS for possible time shifts to determine if one exists
    %     Result  = Conv_X * b_spline_larss_result_2nd_deriv;
    %     %figure;plot(Sim_Ct_larss_Regul_noise,'b');hold on;plot(Result,'r');hold off;
    %     Rms_err = sqrt(sum((Result - Sim_Ct_larss_Regul_noise).^2));
    
    
    % Upsmaple -> shift -> Downsample
    UpSampFactor                     = double(round(min_interval / Upsampling_resolution)) ;
    Sim_AIF_with_noise_Regul_up_samp = interp(Sim_AIF_with_noise_Regul,UpSampFactor);
    
    % Shift the up sampled AIF in wanted times
    time_res_sec     = Upsampling_resolution * 60;
    shift_times      = ( (Min_Time_Delay:time_res_sec:Max_Time_Delay) / 60);
    shift_indices    = round(shift_times/Upsampling_resolution);
    num_shifts       = length(shift_times);
    
    % Initiate shifts matrices results
    CTC_size         = length(Final_Filter_Estimation_Larss);
    Rms_errors_CTC   = zeros(1,num_shifts);
    Est_CTCs         = zeros(CTC_size,num_shifts);
    Rms_errors_BiExp = zeros(1,num_shifts);
    Spline_est       = zeros(CTC_size,num_shifts);
    exp_fit          = zeros(CTC_size,num_shifts);
    
    for i = 1:num_shifts
        
        % Shift for current iteration
        shift_index = shift_indices(i);
        
        %% Create a shifted AIF version (Upsmpl->shift->Downsmp)
        Sim_AIF_with_noise_Regul_up_samp_shifted = zeros(size(Sim_AIF_with_noise_Regul_up_samp));
        
        % Right shift the AIF
        if (shift_index > 0)
            
            Sim_AIF_with_noise_Regul_up_samp_shifted(shift_index+1:end)     = Sim_AIF_with_noise_Regul_up_samp(1:end-shift_index);
            % Duplicate the initial values?
            Sim_AIF_with_noise_Regul_up_samp_shifted(1:shift_index) = Sim_AIF_with_noise_Regul_up_samp_shifted(shift_index);
            
        elseif (shift_index < 0) % Left shift
            
            % Change to a positive index
            shift_index = -1 * shift_index;
            Sim_AIF_with_noise_Regul_up_samp_shifted(1:end-shift_index)    = Sim_AIF_with_noise_Regul_up_samp(shift_index+1:end);
            % Duplicate the last values?
            Sim_AIF_with_noise_Regul_up_samp_shifted(end-shift_index+1:end) = Sim_AIF_with_noise_Regul_up_samp_shifted(end-shift_index);
        else % no shift
            Sim_AIF_with_noise_Regul_up_samp_shifted = Sim_AIF_with_noise_Regul_up_samp;
        end
        % Downsample to original resolution
        Sim_AIF_with_noise_Regul_shifted = downsample(Sim_AIF_with_noise_Regul_up_samp_shifted,UpSampFactor);
        
        %% Create the new Convolution matrix
        [ Conv_X_shift ] = Convolution_Matrix( min_interval, Sim_AIF_with_noise_Regul_shifted );
        
        %% Deconvolve with shifted AIF convolution matrix
        
        if ~Use_Upsampling_Delay_Comp
            
            %             %%
            %             Sim_Ct_larss_Regul_noise_modified = [zeros(length(Sim_Ct_larss_Regul_noise) - 1, 1) ; Sim_Ct_larss_Regul_noise];
            %             time_vec_minutes_modified         = unique([ -1*fliplr(time_vec_minutes) time_vec_minutes]);
            %             knots_modified                    = unique([ -1*fliplr(knots) knots]);
            %             B_mat_modified                    = Create_B_matrix(knots_modified,time_vec_minutes,poly_deg-1);
            %             Conv_X_shift_modified
            
            %%
            B_mat_modified                    = B_mat;
            Conv_X_shift_modified             = Conv_X_shift;
            time_vec_minutes_modified         = time_vec_minutes;
            Sim_Ct_larss_Regul_noise_modified = Sim_Ct_larss_Regul_noise;
            
            
            if RealData_Flag
                [ ~, ~, ~, Spline_est(:,i), ~, ~, ~, idx_fig ] =  ...
                    Regularization_Methods_Simulation( Sim_Ct_larss_Regul, Sim_Ct_larss_Regul_noise_modified, Conv_X_shift_modified, Conv_X_no_noise, time_vec_minutes_modified...
                    , lambda_vec_larss, normalize, min_interval, B_mat_modified, B_PCA, plot_L_Curve, idx_fig, 'Larss', Derivative_Time_Devision, 0, RealData_Flag );
            else
                
                % Deconvolution by regularization for larsson's filter
                [ ridge_regression_result, b_spline_result, b_spline_result_1st_deriv, b_spline_result_2nd_deriv, b_PCA_larss_result, b_PCA_result_1st_deriv, b_PCA_result_2nd_deriv, idx_fig ]...
                    = Regularization_Methods_Simulation(Sim_Ct_larss_Regul, Sim_Ct_larss_Regul_noise_modified,Conv_X_shift_modified,Conv_X_no_noise,time_vec_minutes_modified,...
                    lambda_vec_larss, normalize, min_interval, B_mat_modified, B_PCA, plot_L_Curve, idx_fig , 'Larss' , Derivative_Time_Devision, 0, RealData_Flag);
                
                [ Spline_est(:,i) ] = Choose_Needed_Ht( Filter_Est_Chosen, Sim_Struct.est_larss_filter_Wiener_noise, ridge_regression_result, b_spline_result, b_spline_result_1st_deriv, b_spline_result_2nd_deriv, b_PCA_larss_result, b_PCA_result_1st_deriv, b_PCA_result_2nd_deriv);
            end
            
            %% Estimate Patlak parameters for initial guess
            est_F_noise         = max(Spline_est(:,i));
            if est_F_noise<=0
                % Something is wrong with data. Return the input and stop iterating
                est_delay_by_AIF_correct         = 0;
                Sim_AIF_with_noise_Regul_shifted = Sim_AIF_with_noise_Regul;
                %Final_Filter_Estimation_Larss
                return;
            else
                [est_Ktrans_Patlak_noise, est_Vb_Patlak_noise ,est_E_Patlak_noise, est_MTT_Patlak_noise, idx_fig] = Patlak_Estimation(Sim_Struct, Sim_AIF_with_noise_Regul_shifted, Sim_Ct_larss_Regul_noise_modified,est_F_noise, Verbosity, iter_num, avg_num, idx_fig);
            end
            
            %% Non-linear bi-exp fit
            
            % Prepare parameters for bi-exp fit
            Init_Guess_Larsson    = double( [est_Vb_Patlak_noise est_E_Patlak_noise init_Ve_guess] );
            
            
            if (adjusted_larsson)
                Larsson_function      = @(x,t) Adjusted_Larsson_Filter( t, est_F_noise, x(1), x(2), x(3));
            else
                Larsson_function      = @(x,t) Larsson_Filter( t, est_F_noise, x(1), x(2), x(3), Hct(iter_num));
            end
            
            % Try to avid 0's in data to fit (ydata in optimization problem)
            % tmp          = Spline_est(:,i)/est_F_noise;
            % min_not_zero = min(tmp(tmp>0));
            % toFit        = max(Spline_est(:,i)/est_F_noise, min_not_zero);
            
            if strcmp(FMS_Algorithm,'trust-region-reflective')
                
                % Estimate bi-exp fit to spline result
                [est_params_Larsson_noise,residue_norm_Larsson_noise,residual_Larsson_noise,exitflag_Larsson_noise,algo_info_Larsson_noise] = ...
                    lsqcurvefit(Larsson_function,Init_Guess_Larsson,time_vec_minutes', Spline_est(:,i)/est_F_noise,...
                    LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
                
            elseif strcmp(FMS_Algorithm,'levenberg-marquardt')
                
                Unbounded_Larsson_function = @(x,t) Larsson_function(BoundFunc(x,LowerBound_Larsson,UpperBound_Larsson),t);
                
                [est_params_Larsson_noise,residue_norm_Larsson_noise,residual_Larsson_noise,exitflag_Larsson_noise,algo_info_Larsson_noise] = ...
                    lsqcurvefit(Unbounded_Larsson_function,Init_Guess_Larsson,time_vec_minutes',Spline_est(:,i)/est_F_noise,...
                    [],[],algorithm_options);
                
                est_params_Larsson_noise = BoundFunc(est_params_Larsson_noise,LowerBound_Larsson,UpperBound_Larsson);
                
            elseif strcmp(FMS_Algorithm,'fminsearch')
                
                
                Unbounded_Larsson_function = @(x,t) Larsson_function(BoundFunc(x,LowerBound_Larsson,UpperBound_Larsson),t);
                ydata    = Spline_est(:,i)/est_F_noise;
                RMSCost  = @(vec) sum(vec.^2);
                CostFunc = @(x) RMSCost( Unbounded_Larsson_function(x,time_vec_minutes') - ydata );
                best_transformedParams    = fminsearch(CostFunc, Init_Guess_Larsson, Sim_Struct.algorithm_options);
                est_params_Larsson_noise = BoundFunc(best_transformedParams,LowerBound_Larsson,UpperBound_Larsson);
                
                
                %                 % Estimate bi-exp fit to spline result
                %                 [est_params_Larsson_noise,residue_norm_Larsson_noise,residual_Larsson_noise,exitflag_Larsson_noise,algo_info_Larsson_noise] = ...
                %                     lsqcurvefit(Larsson_function,Init_Guess_Larsson,time_vec_minutes',Spline_est(:,i)/est_F_noise,...
                %                     [],[],algorithm_options);
                
            end
            
            %figure;plot(Sim_Ct_larss_Regul_noise,'b');hold on;plot(Result,'r');hold off;
            %figure;plot(b_spline_larss_result_2nd_deriv);title(['Time Shift: ' num2str(shift_times(i)*60) ]);
            
            Rms_errors_BiExp(i)   = residue_norm_Larsson_noise;
            
            if (adjusted_larsson)
                exp_fit(:,i)          = Adjusted_Larsson_Filter(time_vec_minutes',est_F_noise,est_params_Larsson_noise(1),est_params_Larsson_noise(2),est_params_Larsson_noise(3));
            else
                exp_fit(:,i)          = Larsson_Filter(time_vec_minutes',est_F_noise,est_params_Larsson_noise(1),est_params_Larsson_noise(2),est_params_Larsson_noise(3),Hct(iter_num));
            end
            
            
            
            % Check min RMS for possible time shifts to determine if one exists
            %Result            = Conv_X_shift * b_spline_larss_result_2nd_deriv;
            %Result            = Conv_X * Spline_est(:,i);
            %Result            = Conv_X_shift * Spline_est(:,i);
            Result            = Conv_X_shift * est_F_noise * exp_fit(:,i);
            Est_CTCs(:,i)     = Result;
            
            if RMS_Smooth_Around_Bolus
                
                % Find Bolus Peak
                [peak_val, peak_idx]        = findpeaks(smooth(Sim_Ct_larss_Regul_noise));
                [~, highest_from_peaks_idx] = max(peak_val);
                peak_idx                    = peak_idx(highest_from_peaks_idx);
                %[~, peak_idx]       = max(diff(Sim_Ct_larss_Regul_noise(:,iter_num,avg_num)));
                
                diff_from_bolus_min = Diff_From_Bolus/60; % The difference in seconds(minutes) from the bolus to look on
                diff_from_bolus_idx = round(diff_from_bolus_min/min_interval);
                idx_to_calc         = max(1,peak_idx-diff_from_bolus_idx): (peak_idx+diff_from_bolus_idx);
                
                %                 Sim_Ct_larss_Regul_noise_modified_smooth = smooth(Sim_Ct_larss_Regul_noise_modified);
                %                 A=sum((Est_CTCs(idx_to_calc,:)-repmat(Sim_Ct_larss_Regul_noise_modified_smooth(idx_to_calc),[1 size(Est_CTCs,2)]).^2),1);
                
                Rms_errors_CTC(i)   = sqrt(sum((Result(idx_to_calc) - smooth(Sim_Ct_larss_Regul_noise(idx_to_calc))).^2));
                %Rms_errors_CTC(i)   = sqrt(sum((Result(idx_to_calc) -        Sim_Ct_larss_Regul_noise(idx_to_calc)).^2));
            elseif RMS_Smooth
                Rms_errors_CTC(i) = sqrt(sum((Result - smooth(Sim_Ct_larss_Regul_noise)).^2));
            else
                Rms_errors_CTC(i) = sqrt(sum((Result - Sim_Ct_larss_Regul_noise).^2));
            end
            
            
        else % Up-sampled resolution
            
            UpSampFactor                             = round(min_interval / Upsampling_resolution) ;
            Sim_Ct_larss_Regul_noise_upsmp           = interp(Sim_Ct_larss_Regul_noise,UpSampFactor);
            num_time_stamps                          = length(Sim_AIF_with_noise_Regul_up_samp_shifted);
            
            % Create convolution matrix
            [ Conv_X_shift_upsmp ] = Convolution_Matrix( Upsampling_resolution, Sim_AIF_with_noise_Regul_up_samp_shifted );
            
            % Time vector for AIF and Ct(t)
            time_vec_minutes_upsmp    = (0:num_time_stamps - 1).* Upsampling_resolution;
            %
            B_mat_upsmp = Create_B_matrix(knots,time_vec_minutes_upsmp,poly_deg-1);
            
            
            if RealData_Flag
                [ ~, ~, ~, Spline_est(:,i), ~, ~, ~, idx_fig ] =  ...
                    Regularization_Methods_Simulation( Sim_Ct_larss_Regul, Sim_Ct_larss_Regul_noise_upsmp, Conv_X_shift_upsmp, Conv_X_no_noise, time_vec_minutes_upsmp...
                    , lambda_vec_larss, normalize, min_interval, B_mat_upsmp, B_PCA, plot_L_Curve, idx_fig, 'Larss', Derivative_Time_Devision, 0, RealData_Flag );
            else
                
                % Deconvolution by regularization for larsson's filter            % Deconvolution by regularization for larsson's filter
                [ ridge_regression_result, b_spline_result, b_spline_result_1st_deriv, b_spline_result_2nd_deriv, b_PCA_result_1st_deriv, b_PCA_result_2nd_deriv, idx_fig ]...
                    = Regularization_Methods_Simulation(Sim_Ct_larss_Regul, Sim_Ct_larss_Regul_noise_upsmp,Conv_X_shift_upsmp,Conv_X_no_noise,time_vec_minutes_upsmp,...
                    lambda_vec_larss, normalize, min_interval, B_mat_upsmp, B_PCA, plot_L_Curve, idx_fig , 'Larss' , Derivative_Time_Devision, 0 , RealData_Flag);
                
                [ Spline_est(:,i) ] = Choose_Needed_Ht( Filter_Est_Chosen, est_larss_filter_Wiener_noise, ridge_regression_result, b_spline_result, b_spline_result_1st_deriv, b_spline_result_2nd_deriv, b_PCA_larss_result, b_PCA_result_1st_deriv, b_PCA_result_2nd_deriv);
                
            end
            
            
            % Check min RMS for possible time shifts to determine if one exists
            Result = Conv_X_shift_upsmp * Final_Filter_Estimation_Larss;
            
            %figure;plot(Sim_Ct_larss_Regul_noise,'b');hold on;plot(Result,'r');hold off;
            %figure;plot(b_spline_larss_result_2nd_deriv);title(['Time Shift: ' num2str(shift_times(i)*60) ]);
            
            if  RMS_Smooth_Around_Bolus
                
                % Find Bolus Peak
                [peak_val, peak_idx]        = findpeaks(smooth(Sim_Ct_larss_Regul_noise_upsmp(:,iter_num,avg_num)));
                [~, highest_from_peaks_idx] = max(peak_val);
                peak_idx                    = peak_idx(highest_from_peaks_idx);
                %[~, peak_idx]       = max(diff(Sim_Ct_larss_Regul_noise(:,iter_num,avg_num)));
                
                
                diff_from_bolus_min = Diff_From_Bolus/60; % The difference in seconds(minutes) from the bolus to look on
                diff_from_bolus_idx = round(diff_from_bolus_min/min_interval);
                idx_to_calc         = max(1,peak_idx-diff_from_bolus_idx): (peak_idx+diff_from_bolus_idx);
                
                Rms_errors_CTC(i)   = sqrt(sum((Result(idx_to_calc) - smooth(Sim_Ct_larss_Regul_noise_upsmp(idx_to_calc))).^2));
                %Rms_errors_CTC(i)   = sqrt(sum((Result(idx_to_calc) -        Sim_Ct_larss_Regul_noise_upsmp(idx_to_calc)).^2));
            elseif RMS_Smooth
                Rms_errors_CTC(i) = sqrt(sum((Result - smooth(Sim_Ct_larss_Regul_noise_upsmp)).^2));
            else
                Rms_errors_CTC(i) = sqrt(sum((Result - Sim_Ct_larss_Regul_noise_upsmp).^2));
            end
            
        end
        
    end % For i=1:num_shifts
    
end  % If/Else Simple_AIF_Delay_Correct

% % Plot the erros as a function of time shift
% figure; plot(shift_times*60,Rms_errors_BiExp);
% xlabel('Time shift [sec]');
% ylabel('RMS');
%
% % Check what is the best delay
% [~, min_idx]   = min(Rms_errors_BiExp);
% min_time_shift = shift_times(min_idx);
%
% title(['RMS vs Time Shift. Orig: ' num2str(additional_AIF_delay_sec) ' . Est: ' num2str(min_time_shift*60) ' . Est reg: ' num2str(est_delay_by_spline_result*60) ]);
% figure;
% plot3(shift_times*60,Rms_errors_BiExp,Rms_errors_CTC);
% xlabel('Time Shift [sec]');ylabel('RMS - BiExp');zlabel('RMS - CTC');

% Normalize RMS values
norm_BiExpRms = ( Rms_errors_BiExp - min(Rms_errors_BiExp) ) / max(Rms_errors_BiExp - min(Rms_errors_BiExp)) ;
norm_CTCRms   = ( Rms_errors_CTC - min(Rms_errors_CTC)     ) / max(Rms_errors_CTC   - min(Rms_errors_CTC)  )   ;

% Get the index of the real shift needed
[~, orig_idx] = min(abs(additional_AIF_delay_sec - shift_times*60));

% Find min in both RMS dimensions (average)
%average            = (norm_BiExpRms+norm_CTCRms) / 2;
average            = BiExp2CTC_RMS_Ratio*norm_BiExpRms + (1-BiExp2CTC_RMS_Ratio)*norm_CTCRms;

% Get the first one to be less than 1% error
Zero_shift_index          = round(length(average)/2);
indices                   = find(average < 0.01);
[~ , first_min_error_idx] = min(abs(indices-Zero_shift_index));
min_idx                   = indices(first_min_error_idx);

% [~, min_idx]       = min(average);
% min_idx            = min_idx(1); % Take the first one if there is more than 1

[~, min_idx_CTC]   = min(norm_CTCRms);
min_idx_CTC        = min_idx_CTC(1); % Take the first one if there is more than 1

[~, min_idx_BiExp] = min(norm_BiExpRms);
min_idx_BiExp      = min_idx_BiExp(1); % Take the first one if there is more than 1

% Estimate the delay
min_shift_sec_BiExp           = shift_times(min_idx_BiExp)*60;
min_shift_sec_CTC             = shift_times(min_idx_CTC)*60;
est_delay_by_AIF_correct      = shift_times(min_idx)*60;


% CTC = smooth(Sim_Ct_larss_Regul_noise);
% CTC = Sim_Ct_larss_Regul_noise;
% figure;
% plot(Est_CTCs(:,min_idx_CTC),'r');
% hold on;
% plot(Est_CTCs(:,orig_idx),'g');
% plot(CTC,'Color','k','LineWidth',5);
% hold off;


if (plot_flag)
    
    fig_num = figure;
    
    subplot(3,1,1);
    plot(shift_times*60,norm_BiExpRms);
    xlabel('Time Shift [sec]');
    ylabel('RMS - BiExp');
    title(['Orig Time Shift : ' num2str(additional_AIF_delay_sec) ' .Est. - weighted RMSs BiExp : ' num2str(min_shift_sec_BiExp) ' [Sec]']);
    
    subplot(3,1,2);
    plot(shift_times*60,norm_CTCRms);
    xlabel('Time Shift [sec]');
    ylabel('RMS - CTC');
    title(['Orig Time Shift : ' num2str(additional_AIF_delay_sec) ' .Est. - weighted RMSs CTC : ' num2str(min_shift_sec_CTC) ' [Sec]']);
    
    
    % plot3(shift_times*60,norm_BiExpRms,norm_CTCRms);
    % xlabel('Time Shift [sec]');ylabel('RMS - BiExp');zlabel('RMS - CTC');
    % title(['Orig Time Shift : ' num2str(additional_AIF_delay_sec) ' .Est. - weighted RMSs : ' num2str(min_shift_sec) ' [Sec]']);
    % xlabel('Time Shift [sec]');
    % ylabel('Bi-Exp fit Error');
    % % Display the chosen point
    % hold on;
    % plot3(shift_times(min_idx)*60,norm_BiExpRms(min_idx),norm_CTCRms(min_idx),'m*','MarkerSize',20);
    % hold off;
    
    
    subplot(3,1,3);
    h1 = plot(time_vec_minutes,Sim_Ct_larss_Regul_noise,'b*');
    hold on;
    if RMS_Smooth_Around_Bolus
        plot(time_vec_minutes(idx_to_calc),Sim_Ct_larss_Regul_noise(idx_to_calc),'co');
    end
    smooth_CTC = smooth(Sim_Ct_larss_Regul_noise);
    h2 = plot(time_vec_minutes,smooth_CTC,'--m');
    %plot(time_vec_minutes(idx_to_calc),smooth_CTC(idx_to_calc),'co');
    
    h3 = plot(time_vec_minutes, Est_CTCs(:,orig_idx),'gx');
    h4 = plot(time_vec_minutes, Est_CTCs(:,min_idx),'ro');
    h5 = plot(time_vec_minutes, Est_CTCs(:,min_idx_CTC),'kx');
    %h3 = plot(time_vec_minutes, Conv_X_shift * Spline_est(:,orig_idx),'gx');
    %h4 = plot(time_vec_minutes, Conv_X_shift * Spline_est(:,min_idx),'ro');
    %h4 = plot(time_vec_minutes, Conv_X_shift * Spline_est(:,min_idx_CTC),'ro');
    hold off;
    xlabel('Time [Min]');
    ylabel('Amplitude');
    legend([h1 h2 h3 h4 h5], 'CTC', 'Smooth CTC', 'Fitted - Orig. Delay', 'Fitted - Est. Weighted', 'Fitted - Est. CTC RMS');
    %figure;plot(b_spline_larss_result_2nd_deriv);title(['Time Shift: ' num2str(shift_times(i)*60) ]);
    title('Orig. CTC, Fitted splines (est. and orig. delay)');
    
    % subplot(3,1,3);
    % title('Bi-Exp fit - Original h(t) vs. Estimated');
    
    %     % Printing image to PDF
    %     gprint(fig_num,'Run_Output/Time_Shift_Estimation.png');
    %
    %     idx_fig    = idx_fig + 1;
    %     idx_string = ['idx_' num2str(idx_fig,'%03i')];
    %     AddToLog('Run_Output/',idx_string,'\\subsection*{\\underline{AIF Delay Estimation}}');
    %     idx_fig    = idx_fig + 1;
    %     idx_string = ['idx_' num2str(idx_fig,'%03i')];
    %     AddToLog('Run_Output/',idx_string,'AIFDelayEst','Time_Shift_Estimation.png');
    %
    % Print result to PDF
    [Sim_Struct.idx_fig] = Print2Pdf(fig_num, Sim_Struct.idx_fig, 'Time_Shift_Estimation.png', './Run_Output/', 'AIF Delay Estimation', 'AIFDelayEst');
    
end

%display('');figure;plot(time_vec_minutes,Spline_est(:,28));hold on; plot(time_vec_minutes,exp_fit(:,28),'r');hold off;

% Assign the new spline after fixing for delay
Final_Filter_Estimation_Larss = Spline_est(:,min_idx);


end

