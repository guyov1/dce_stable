function [ t_delay, sigma, Est_ht, conv_result_ht, RMS_ht, RMS_gauss, RMS_params ] = Get_Ht_Deconvolving( AIF, Ct , sec_interval, DEBUG)
%Get_Ht_Deconvolving Extracting possible h_t filter
%   The function gets AIF(nT), Ct(nT) - the arterial input function and the
%   tissue conventration (in minutes metric).
%   In addition it gets the time interval between samples (sec_interval).
%   It uses deconvolution (Wiener Filter) to estimate the convoultion
%   filter ( h(t) ).
%   The SNR approximation is based on experimental results.
%   After h(t) aprox. , it fits a gaussian using non-linear fitting method.
%   The output consist of:
%   - t_delay and sigma, the gaussian parameters
%   - h_t, the deconvolved filter estimation
%   - conv_result, the result of AIF convolved with h_t
%   - RMS_ht, the root mean square error of conv_result with h_t and Ct.
%   - RMS_gauss, the root mean square error of conv_result with estimated gaussian and Ct.
%   - RMS_params, the root mean square error of the gaussian using the parameters
%          and h(t)  (if too big, h(t) is not really a gaussian)

%% Initiate parameters

% Time interval between samples
Fs             = 1 / sec_interval; %[Hz] - Sampling rate
min_interval   = sec_interval/60;  %[min]

% Get number of time stamps (vector length)
if ( size(Ct,1) > 1)
    num_time_stamps    = size(Ct,1);
else
    num_time_stamps    = size(Ct,2);
end

% Time vector for AIF and Ct(t)
time_vec_minutes   = (1:num_time_stamps).* min_interval;

% Non-linear square fitting parameters
FMS_TolFun      = 1e-11;
FMS_MaxFunEvals = 10000;
FMS_MaxIter     = 10000;
options         = optimset('TolFun',FMS_TolFun,'MaxFunEvals',FMS_MaxFunEvals,'MaxIter',FMS_MaxIter,'Display','off');

% Initial Guess (td=1sec, var=1)
init_guess = [1 1];

% Set lower and upper bounds for parameters
% Time delay -> 0 to 3 seconds
% Sigma      -> 0 to 1 seconds
LowerBound  = [0 0];
UpperBound  = [3 1];

%% Estimating h(t) by Wiener filter
[Est_ht] = Wiener_Filter( AIF, Ct, Fs);
Est_ht_T = transpose(Est_ht);

% The analytic funcational of a gaussian function
Gaussian_function = @(x,t) Gaussian( t,x(1),x(2) );

% lsqcurvefit parameters are:
% analytic function, initial parameters, time vector, data points ,lower
% and upper bounds and algorithm options
[est_params,residue_norm,residual,exitflag,algo_info] = ...
    lsqcurvefit(Gaussian_function,init_guess,time_vec_minutes',Est_ht_T,LowerBound,UpperBound,options);

% Put parameters in wanted variables
t_delay_in_min = est_params(1);
t_delay        = 60 * t_delay_in_min; % Convert delay time to seconds
est_var        = est_params(2);
sigma_in_min   = sqrt(est_var);
sigma          = 60 * sigma_in_min; % Convert sigma to seconds

% Calculate gaussian outo of parameters
calculated_gaussian = Gaussian(time_vec_minutes, t_delay_in_min, sigma_in_min^2);

% Calculate RMS of convolution result comparing to Ct(t)
RMS_params = sqrt( sum( (calculated_gaussian - Est_ht).^2 ) );

%% Filter AIF through kernel

% Filter the AIF with the estimated ht
conv_result_ht       = filter(Est_ht,1,AIF);
% Filter the AIF with the gaussian kernel
conv_result_gaussian = filter(calculated_gaussian,1,AIF);

% Zero negative values
conv_result_ht(conv_result_ht<0)             = 0;
conv_result_gaussian(conv_result_gaussian<0) = 0;

% Calculate RMS of convolution results comparing to Ct(t)
RMS_ht    = sqrt( sum( (Ct - conv_result_ht).^2 ) );
RMS_gauss = sqrt( sum( (Ct - conv_result_gaussian).^2 ) );

%% Debugging results

if (DEBUG)
    
    figure;
    
    % Plot estimated h(t) and calculated gaussian
    %hold on;plot(Est_ht,'b');plot(calculated_gaussian,'g');hold off;
    subplot(1,3,1);
    plot(time_vec_minutes,AIF,time_vec_minutes,AIF,'*');
    title('Input AIF');
    xlabel('Time [Min]');
    legend('Input AIF');
    
    subplot(1,3,2);
    hold on;
    h1 = plot(time_vec_minutes,Est_ht,'b');
    h2 = plot(time_vec_minutes,Est_ht,'g*');
    h3 = plot(time_vec_minutes,calculated_gaussian,'r');
    h4 = plot(time_vec_minutes,calculated_gaussian,'y*');
    hold off;
    title('Estimated Filter h(t)');
    xlabel('Time [Min]')
    legend([h2 h4],'Estimated h(t)','Fitted Gaussian');
    
    
    subplot(1,3,3);
    hold on;
    h1 = plot(time_vec_minutes,Ct,'b');
    h2 = plot(time_vec_minutes,Ct,'g*');
    h3 = plot(time_vec_minutes,conv_result_ht,'r');
    h4 = plot(time_vec_minutes,conv_result_ht,'y*');
    h5 = plot(time_vec_minutes,conv_result_gaussian,'k');
    h6 = plot(time_vec_minutes,conv_result_gaussian,'m*');
    hold off;
    title('Input Ct and est. Ct');
    xlabel('Time [Min]');
    legend([h2 h4 h6],'Input Ct(t)','Estimated Ct(t) - h(t)','Estimated Ct(t) - Gaussian');
    
end

end

