function [ est_delay_by_AIF_correct, Sim_AIF_with_noise_Regul_shifted,b_spline_larss_result_2nd_deriv, idx_fig ] = AIF_Delay_Correct( Sim_Struct, ht_Struct, Verbosity, iter_num, avg_num, idx_fig)

% Take from struct variables used in local function
normalize                       = Sim_Struct.normalize;
B_mat                           = Sim_Struct.B_mat;
B_PCA                           = Sim_Struct.PCA_B_mat;
plot_L_Curve                    = Sim_Struct.plot_L_Curve;
Derivative_Time_Devision        = Sim_Struct.Derivative_Time_Devision;
lambda_vec_larss                = Sim_Struct.lambda_vec_larss;
min_interval                    = Sim_Struct.min_interval(iter_num);
time_vec_minutes                = Sim_Struct.time_vec_minutes;
Upsampling_resolution           = Sim_Struct.Upsampling_resolution;
Max_Time_Delay                  = Sim_Struct.Max_Time_Delay;
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

b_spline_larss_result_2nd_deriv = ht_Struct.b_spline_larss_result_2nd_deriv;
Sim_AIF_with_noise_Regul        = ht_Struct.Sim_AIF_with_noise_Regul;
Sim_Ct_larss_Regul              = ht_Struct.Sim_Ct_larss_Regul;
Sim_Ct_larss_Regul_noise        = ht_Struct.Sim_Ct_larss_Regul_noise;
Conv_X_no_noise                 = ht_Struct.Conv_X_no_noise;

% ---------------------------------------------------------------------
%                 Regular correction by looking on h(t) shift
% ---------------------------------------------------------------------
est_delay_by_spline_result       = NaN; % Initiate with NaN
if Sim_Struct.Simple_AIF_Delay_Correct
    
    [~, max_idx] = max(b_spline_larss_result_2nd_deriv);
    
    % If the maximum is not the first value, we might think there is a delay
    if (max_idx ~= 1)
        
        [peak_val, peak_idx]       = findpeaks(b_spline_larss_result_2nd_deriv);
        
        % We expect a peak in case of delay
        if ~isempty(peak_idx)
            
            % The estimated delay in minutes
            est_delay_by_spline_result                   = time_vec_minutes(peak_idx(1));
            
            % Shift the AIF according to estimation
            shift_index                                             = peak_idx(1);
            Sim_AIF_with_noise_Regul_shifted                        = zeros(size(Sim_AIF_with_noise_Regul));
            Sim_AIF_with_noise_Regul_shifted(1:end-shift_index)     = Sim_AIF_with_noise_Regul(shift_index+1:end);
            Sim_AIF_with_noise_Regul_shifted(end-shift_index+1:end) = Sim_AIF_with_noise_Regul_shifted(end-shift_index); % Duplicate the last values?
            
            % Create new convolution matrix for
            [ Conv_X_shifted ] = Convolution_Matrix( min_interval, Sim_AIF_with_noise_Regul_shifted );
            
            [ridge_regression_larss_result, b_spline_larss_result, b_spline_larss_result_1st_deriv, b_spline_larss_result_2nd_deriv, b_PCA_gauss_result_2nd_deriv, idx_fig]...
                = Regularization_Methods_Simulation(Sim_Ct_larss_Regul, Sim_Ct_larss_Regul_noise,Conv_X_shifted,Conv_X_no_noise,time_vec_minutes,...
                lambda_vec_larss, normalize, min_interval, B_mat, plot_L_Curve, idx_fig , 'Larss' , Derivative_Time_Devision, plot_flag );
            
        end
        
    end
    
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
    UpSampFactor                     = round(min_interval / Upsampling_resolution) ;
    Sim_AIF_with_noise_Regul_up_samp = interp(Sim_AIF_with_noise_Regul,UpSampFactor);
    
    % Shift the up sampled AIF in wanted times
    time_res_sec     = Upsampling_resolution * 60;
    shift_times      = ( (-Max_Time_Delay:time_res_sec:Max_Time_Delay) / 60);
    shift_indices    = round(shift_times/Upsampling_resolution);
    num_shifts       = length(shift_times);
    
    % Initiate shifts matrices results
    CTC_size         = length(b_spline_larss_result_2nd_deriv);
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
            Sim_Ct_larss_Regul_modified       = Sim_Ct_larss_Regul;
            Conv_X_shift_modified             = Conv_X_shift;
            time_vec_minutes_modified         = time_vec_minutes;
            Sim_Ct_larss_Regul_noise_modified = Sim_Ct_larss_Regul_noise;
            
        
            % Deconvolution by regularization for larsson's filter
            [~, ~, ~, ~, Spline_est(:,i), idx_fig]...
                = Regularization_Methods_Simulation(Sim_Ct_larss_Regul, Sim_Ct_larss_Regul_noise_modified,Conv_X_shift_modified,Conv_X_no_noise,time_vec_minutes_modified,...
                lambda_vec_larss, normalize, min_interval, B_mat_modified, B_PCA, plot_L_Curve, idx_fig , 'Larss' , Derivative_Time_Devision, 0 );
            
            %% Estimate Patlak parameters for initial guess
            AIF                   = Sim_AIF_with_noise_Regul_shifted;
            Y_vec_Vb              = Sim_Ct_larss_Regul ./ AIF;                        %[mL/100g]
            X_vec                 = cumtrapz(time_vec_minutes,AIF) ./ AIF; %[min]
            [~, bolus_idx]        = max(abs(diff(AIF)));
            base_value            = mean(AIF(1:max(bolus_idx-1,1)));
            mult_val_Thresh       = 3;
            Threshold             = mult_val_Thresh*base_value;
            stable_idx            = find(AIF > Threshold);
            % If threshold was too big, decrease it till we get some value
            while (isempty(stable_idx)|| length(stable_idx)<2)
                mult_val_Thresh       = mult_val_Thresh/1.5;
                Threshold             = mult_val_Thresh*base_value;
                stable_idx            = find(AIF > Threshold);
            end
            % Take the stable points out of the vector
            X_vec                 = X_vec(stable_idx);
            Y_vec_Vb              = Y_vec_Vb(stable_idx);
            % Remove Zeros/NaNs/Infs caused because of division by 0
            nan_indices           = find(isnan(Y_vec_Vb));
            inf_indices           = find(~isfinite(Y_vec_Vb));
            Y_vec_Vb(nan_indices) = [];
            Y_vec_Vb(inf_indices) = [];
            X_vec(nan_indices)    = [];
            X_vec(inf_indices)    = [];
            % Fine straight line coefficent
            [linear_params]       = polyfit(X_vec,Y_vec_Vb,1);
            % Ki,Vb estimation
            est_Ki_Patlak_noise   = linear_params(1); %a in ax+b
            est_Vb_Patlak_noise   = linear_params(2); %b in ax+b
            
            %% Non-linear bi-exp fit
            
            % Prepare parameters for bi-exp fit
            estF                  = max(Spline_est(:,i));
            E_Patlak_est          = est_Ki_Patlak_noise / estF; % E = Ki / F
            Init_Guess_Larsson    = double( [est_Vb_Patlak_noise E_Patlak_est 10] );
            Larsson_function      = @(x,t) Larsson_Filter( t, estF, x(1), x(2), x(3), Hct(iter_num));
            
            % Estimate bi-exp fit to spline result
            [est_params_Larsson_noise,residue_norm_Larsson_noise,residual_Larsson_noise,exitflag_Larsson_noise,algo_info_Larsson_noise] = ...
                lsqcurvefit(Larsson_function,Init_Guess_Larsson,time_vec_minutes',Spline_est(:,i)/estF,...
                LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
            
            %figure;plot(Sim_Ct_larss_Regul_noise,'b');hold on;plot(Result,'r');hold off;
            %figure;plot(b_spline_larss_result_2nd_deriv);title(['Time Shift: ' num2str(shift_times(i)*60) ]);
            
            Rms_errors_BiExp(i)   = residue_norm_Larsson_noise;
            exp_fit(:,i)          = Larsson_Filter(time_vec_minutes',estF,est_params_Larsson_noise(1),est_params_Larsson_noise(2),est_params_Larsson_noise(3),Hct(iter_num));
            
            % Check min RMS for possible time shifts to determine if one exists
            %Result            = Conv_X_shift * b_spline_larss_result_2nd_deriv;
            %Result            = Conv_X * Spline_est(:,i);
            %Result            = Conv_X_shift * Spline_est(:,i);
            Result            = Conv_X_shift * estF * exp_fit(:,i);
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
            
            % Deconvolution by regularization for larsson's filter
            [~, ~, ~, ~, Spline_est(:,i), idx_fig]...
                = Regularization_Methods_Simulation(Sim_Ct_larss_Regul, Sim_Ct_larss_Regul_noise_upsmp,Conv_X_shift_upsmp,Conv_X_no_noise,time_vec_minutes_upsmp,...
                lambda_vec_larss, normalize, min_interval, B_mat_upsmp, plot_L_Curve, idx_fig , 'Larss' , Derivative_Time_Devision, 0 );
            
            % Check min RMS for possible time shifts to determine if one exists
            Result = Conv_X_shift_upsmp * b_spline_larss_result_2nd_deriv;
            
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
b_spline_larss_result_2nd_deriv = Spline_est(:,min_idx);


end

