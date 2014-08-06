function [ b_spline_result_2nd_deriv, est_delay ] = Regularization_Methods( Sim_Ct_T_noise, AIF_noise, Conv_Matrix,time_vec_minutes, lambda_vec, normalize, min_interval, B_mat, Delay_Correct)
% Regularization_Methods - Implement deconvolution by regularization
%
% Formalizing the problem:
%
%   Without splines:
%
%   A*R_vec = C_vec
%   min { ||A*R_vec - C_vec||^2  + lambda^2||L*R_vec ||^2 }
%
%   where -
%            A     - Conv_Matrix
%            R_vec - Model vector ( exp. for Larson - F*IRF(t) )
%            C_vec - Ct(t) - Sim_Ct_T_noise - Tissue concentration
%            L     - Some matrix ( Identity/1st Derivative/2nd Derivative)
%
%   With splines:
%   Force the solution R_vec to be in a splines form. Hence,
%   R_vec = B*V_vec
%
%   where -
%            B     - Spline matrix
%            V_vec - Polynomial coefficients for B splines
%
%   So the new minimzation problem, takes the form:
%
%   A*B*V_vec = C_vec
%   If we define A*B = D, we can rewrite:
%   D*V_vec = C_vec
%   min { ||D*V_vec - C_vec||^2  + lambda^2||L*V_vec ||^2 }

% Using derivative regularization constraint
num_points   = max(size(time_vec_minutes));

% Second derivative
%deriv_matrix   = (1/min_interval)*toeplitz([-1,zeros(1,size(B_mat,2)-1)],[[-1 1],zeros(1,size(B_mat,2)-2)]);
deriv_matrix   = toeplitz([-1,zeros(1,size(B_mat,2)-1)],[[-1 1],zeros(1,size(B_mat,2)-2)]);
deriv_matrix_2 = deriv_matrix*deriv_matrix;

A_new_2        = Conv_Matrix*B_mat*inv(deriv_matrix_2);

if (normalize == 1)
    knots_coeff_3 = ridge(Sim_Ct_T_noise,A_new_2,lambda_vec(4),0);
    knots_coeff_3 = knots_coeff_3(2:end);
else
    knots_coeff_3 = ridge(Sim_Ct_T_noise,A_new_2,lambda_vec(4),1);
end

inv_matrix_2              = B_mat * inv(deriv_matrix_2);
b_spline_result_2nd_deriv = inv_matrix_2 * knots_coeff_3;


% Remove zeros
b_spline_result_2nd_deriv = max(b_spline_result_2nd_deriv,0);


% Check for possible AIF delay
if (Delay_Correct)
    
    % Delay estimation parameters
    Upsampling_resolution   = 0.1 / 60; % Set the upsampling target
    Max_Time_Delay          = 3;
    RMS_Smooth_Around_Bolus = false;
    RMS_Smooth              = true;
    
    % Create convolution indices
    nTimePoints = size(b_spline_result_2nd_deriv,1);
    [X, Y]      = meshgrid(1:nTimePoints);
    Out         = ((Y-X)+1).*(X<=Y);
    Out(Out==0) = nTimePoints+1; % The zero places will hold NaN temporarily
    
    knot_interval = 2;
    knots         = time_vec_minutes(1:knot_interval:end);
    poly_deg      = 4;
    
    plot_L_Curve             = false;
    idx_fig                  = 1;
    Derivative_Time_Devision = false;
    
    lambda_vec_larss = [4.7 0.9 0.1 0.2];
    
    LowerBound_Larsson  = [0   0     0];
    UpperBound_Larsson  = [100 1     100];
    
    Hct                 = 0.38;
    
    Diff_From_Bolus                        = 10;        % The difference in seconds from the bolus to look on
    
    
    FMS_TolFun          = 1e-11;
    FMS_MaxFunEvals     = 10000;
    FMS_MaxIter         = 10000;
    %FMS_Algorithm      = 'levenberg-marquardt'; % Does not supper lower-upper bounds
    FMS_Algorithm       = 'trust-region-reflective';
    algorithm_options   = optimset('TolFun',FMS_TolFun,'MaxFunEvals',FMS_MaxFunEvals,'MaxIter',FMS_MaxIter,...
        'Display','off','Algorithm',FMS_Algorithm);
    
    % Upsmaple -> shift -> Downsample
    UpSampFactor                     = double(round(min_interval / Upsampling_resolution)) ;
    Sim_AIF_with_noise_Regul_up_samp = interp(AIF_noise,UpSampFactor);
    
    Sim_Ct_larss_Regul_noise = Sim_Ct_T_noise;
    
    % Shift the up sampled AIF in wanted times
    time_res_sec     = Upsampling_resolution * 60;
    shift_times      = ( (-Max_Time_Delay:time_res_sec:Max_Time_Delay) / 60);
    shift_indices    = round(shift_times/Upsampling_resolution);
    num_shifts       = length(shift_times);
    
    % Initiate shifts matrices results
    CTC_size         = length(b_spline_result_2nd_deriv);
    Rms_errors_CTC   = zeros(1,num_shifts);
    Est_CTCs         = zeros(CTC_size,num_shifts);
    Rms_errors_BiExp = zeros(1,num_shifts);
    Spline_est       = zeros(CTC_size,num_shifts);
    exp_fit          = zeros(CTC_size,num_shifts);
    
    % Try different possible shifts
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
        
        B_mat_modified                    = B_mat;
        Sim_Ct_larss_Regul_modified       = Sim_Ct_larss_Regul_noise;
        Conv_X_shift_modified             = Conv_X_shift;
        time_vec_minutes_modified         = time_vec_minutes;
        Sim_Ct_larss_Regul_noise_modified = Sim_Ct_larss_Regul_noise;
        
        Sim_Ct_larss_Regul                = Sim_Ct_larss_Regul_noise_modified;
        
        Conv_X_no_noise = Conv_X_shift_modified;
        
        % Deconvolution by regularization for larsson's filter
        [~, ~, ~, Spline_est(:,i), idx_fig]...
            = Regularization_Methods_Simulation(Sim_Ct_larss_Regul_modified, Sim_Ct_larss_Regul_noise_modified,Conv_X_shift_modified,Conv_X_no_noise,time_vec_minutes_modified,...
            lambda_vec_larss, normalize, min_interval, B_mat_modified, plot_L_Curve, idx_fig , 'Larss' , Derivative_Time_Devision, 0 );

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
        Larsson_function      = @(x,t) Larsson_Filter( t, estF, x(1), x(2), x(3), Hct);
        
        % Estimate bi-exp fit to spline result
        [est_params_Larsson_noise,residue_norm_Larsson_noise,residual_Larsson_noise,exitflag_Larsson_noise,algo_info_Larsson_noise] = ...
            lsqcurvefit(Larsson_function,Init_Guess_Larsson,time_vec_minutes',Spline_est(:,i)/estF,...
            LowerBound_Larsson,UpperBound_Larsson,algorithm_options);
        
        %figure;plot(Sim_Ct_larss_Regul_noise,'b');hold on;plot(Result,'r');hold off;
        %figure;plot(b_spline_larss_result_2nd_deriv);title(['Time Shift: ' num2str(shift_times(i)*60) ]);
        
        Rms_errors_BiExp(i)   = residue_norm_Larsson_noise;
        exp_fit(:,i)          = Larsson_Filter(time_vec_minutes',estF,est_params_Larsson_noise(1),est_params_Larsson_noise(2),est_params_Larsson_noise(3),Hct);
        
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
        
        
    end % For i=1:num_shifts
    
    % Normalize RMS values
    norm_CTCRms   = ( Rms_errors_CTC - min(Rms_errors_CTC)     ) / max(Rms_errors_CTC   - min(Rms_errors_CTC)  )   ;
    
    % Find min in both RMS dimensions (average)
    average            = norm_CTCRms;
    
    % Get the first one to be less than 1% error
    Zero_shift_index          = round(length(average)/2);
    indices                   = find(average < 0.01);
    [~ , first_min_error_idx] = min(abs(indices-Zero_shift_index));
    min_idx                   = indices(first_min_error_idx);
    
    % [~, min_idx]       = min(average);
    % min_idx            = min_idx(1); % Take the first one if there is more than 1
    
    [~, min_idx_CTC]   = min(norm_CTCRms);
    min_idx_CTC        = min_idx_CTC(1); % Take the first one if there is more than 1
    est_delay_by_AIF_correct      = shift_times(min_idx)*60;
    
    est_delay = est_delay_by_AIF_correct;
    
    % Assign the new spline after fixing for delay
    b_spline_result_2nd_deriv = Spline_est(:,min_idx);
    
else
    est_delay = 0;
end


end

