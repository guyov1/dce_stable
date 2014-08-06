% Generate a 4 parameters uniformly

% Larsson filter parameters
F_low_val  = 5;
F_max_val  = 140;
Vb_low_val = 3;
Vb_max_val = 20;
Ve_low_val = 1;
Ve_max_val = 10;
E_low_val  = 0;
E_max_val  = 1;

Hct        = 0.38;

Num_iterations = 1000;

F_vec      = F_low_val  + (F_max_val  - F_low_val) * rand(Num_iterations,1);
Vb_vec     = Vb_low_val + (Vb_max_val - Vb_low_val)* rand(Num_iterations,1);
Ve_vec     = Ve_low_val + (Ve_max_val - Vb_low_val)* rand(Num_iterations,1);
E_vec      = E_low_val  + (E_max_val  - E_low_val) * rand(Num_iterations,1);

sec_interval_low_res     = 2;
sec_interval_high_res    = 0.01;
min_interval_low_res     = sec_interval_low_res/60;
min_interval_high_res    = sec_interval_high_res/60;
total_sim_time_min       = 4;
num_time_stamps_high_res = round(total_sim_time_min/min_interval_high_res);
num_time_stamps_low_res  = round(total_sim_time_min/min_interval_low_res);
time_vec_high_resolution = (0:num_time_stamps_high_res-1) * min_interval_high_res;
time_vec_sub_sampled     = (0:num_time_stamps_low_res-1) * min_interval_low_res;
UpSampFactor             = round(min_interval_low_res / min_interval_high_res) ;

% Create AIF
% AIF average population parameters
A1    = 0.809;
A2    = 0.330;
% Changed to shift the boluses to the right
T1    = 0.17046;
T2    = 0.365;
%T1    = 0.27046;
%T2    = 0.465;
sig1  = 0.0563;
sig2  = 0.132;
alpha = 1.050;
beta  = 0.1685;
s     = 38.078;
tau   = 0.483;
Sim_AIF_high_res_no_noise              = AIF_Parker(time_vec_high_resolution,A1,sig1,T1,A2,sig2,T2,alpha,beta,s,tau,0)'; %[mM]
SNR_ratio                              = 15;
noise_sigma_AIF                        = mean(Sim_AIF_high_res_no_noise) ./ SNR_ratio;
noise_to_add_AIF                       = noise_sigma_AIF .* randn(num_time_stamps_high_res,1);
Sim_AIF_high_res_with_noise            = Sim_AIF_high_res_no_noise + noise_to_add_AIF;

% Initiate error vector
J = zeros(Num_iterations,1);

% Initiate parameters estimation
est_F_vec  = zeros(Num_iterations,1);
est_Vb_vec = zeros(Num_iterations,1);
est_E_vec  = zeros(Num_iterations,1);
est_Ve_vec = zeros(Num_iterations,1);

% Create convolution indices
[ Conv_AIF ] = Convolution_Matrix( min_interval(iter_num), Sim_AIF_high_res_with_noise );

% Function to be fitted
CTC_function_4_Sourbron  = @(x,t) CTC_function( t, Conv_AIF, x(1), x(2), x(3), x(4), Hct);
CTC_function_3_Params    = @(x,t) CTC_function( t, Conv_AIF, F_orig, x(2), x(3), x(4), Hct);

% Gues initial F randomly
Vb_init             = Vb_low_val + (Vb_max_val - Vb_low_val)* rand;
Ve_init             = Ve_low_val + (Ve_max_val - Vb_low_val)* rand;
E_init              = E_low_val  + (E_max_val  - E_low_val) * rand;
Init_Guess_Larsson  = [Vb_init E_init Ve_init];
F_init_guess        = F_low_val  + (F_max_val  - F_low_val) * rand;
Init_Guess_Sourbron = [F_init_guess Init_Guess_Larsson];

LowerBound_Larsson  = [0   0     0];
UpperBound_Larsson  = [Vb_max_val E_max_val     Ve_max_val];
LowerBound_Sourbron = [F_low_val LowerBound_Larsson];
UpperBound_Sourbron = [F_max_val UpperBound_Larsson];


% Non-linear square fitting parameters
FMS_TolFun          = 1e-11;
FMS_MaxFunEvals     = 10000;
FMS_MaxIter         = 10000;
%FMS_Algorithm      = 'levenberg-marquardt'; % Does not supper lower-upper bounds
FMS_Algorithm       = 'trust-region-reflective';
algorithm_options   = optimset('TolFun',SimParams.FMS_TolFun,'MaxFunEvals',SimParams.FMS_MaxFunEvals,'MaxIter',SimParams.FMS_MaxIter,...
    'Display','off','Algorithm',SimParams.FMS_Algorithm);

% For each set of parameters, calculate the RMS error
for i = 1 : Num_iterations
    
    Filter_high_res             = F_vec(i)*Larsson_Filter(time_vec_high_resolution, F_vec(i), Vb_vec(i), E_vec(i), Ve_vec(i), Hct);
    Filter_low_res              = downsample(Filter_high_res,UpSampFactor);
    
    Sim_Ct_high_res             = filter(Filter_high_res*min_interval_high_res,1,Sim_AIF_high_res_no_noise);
    noise_sigma_Ct              = mean(Sim_Ct_high_res) ./ SNR_ratio;
    noise_to_add_Ct             = noise_sigma_Ct .* randn(num_time_stamps_high_res,1);
    Sim_Ct_high_res_with_noise  = Sim_Ct_high_res + noise_to_add_Ct;
    
    Sim_Ct_low_res              = downsample(Sim_Ct_high_res,UpSampFactor);
    Sim_Ct_low_res_noise        = downsample(Sim_Ct_high_res_with_noise,UpSampFactor);
    
    % Estimating 4 paramaeters
    % F, Vb, E, Ve
    [est_params_Sourbron_noise,residue_norm,residual,exitflag,algo_info] = ...
        lsqcurvefit(CTC_function_4_Sourbron,Init_Guess_Sourbron,time_vec_sub_sampled',Sim_Ct_low_res_noise,...
        LowerBound_Sourbron,UpperBound_Larsson,algorithm_options);
    % Estimating 3 paramaeters
    % Vb, E, Ve
    [est_params_Sourbron_noise,residue_norm,residual,exitflag,algo_info] = ...
        lsqcurvefit(CTC_function_3_Params,Init_Guess_Larsson,time_vec_sub_sampled',F_vec(i),Sim_Ct_low_res_noise,...
        LowerBound_Larsson,UpperBound_Sourbron,algorithm_options);
    
    if length(est_params_Sourbron_noise) == 4
        % 3 parameters estimation
        est_F_vec(i)  = est_params_Sourbron_noise(1);
        est_Vb_vec(i) = est_params_Sourbron_noise(2);
        est_E_vec(i)  = est_params_Sourbron_noise(3);
        est_Ve_vec(i) = est_params_Sourbron_noise(4);
        
        J(i)   = residue_norm;
        
    elseif length(est_params_Sourbron_noise) == 3
        % 4 parameters estimation
        est_F_vec(i)  = F_vec(i);
        est_Vb_vec(i) = est_params_Sourbron_noise(1);
        est_E_vec(i)  = est_params_Sourbron_noise(2);
        est_Ve_vec(i) = est_params_Sourbron_noise(3);
        
    end  
    
    est_filter = est_F_vec(i)*Larsson_Filter(time_vec_sub_sampled, est_F_vec(i), est_Vb_vec(i), est_E_vec(i), est_Ve_vec(i), Hct);
    
    % Plot the estimated parameteric filter vs. the original
    figure;
    plot(time_vec_sub_sampled,est_filter,'r-');
    hold on;
    plot(time_vec_sub_sampled,Filter_low_res,'g');
    hold off;
    title('Estimated parameteric filter vs. original');
    
    % Fitted CTC
    Fitted_CTC = CTC_function( time_vec_sub_sampled', Conv_AIF, est_F_vec(i), est_Vb_vec(i), est_E_vec(i), est_Ve_vec(i), Hct);
    figure;
    plot(time_vec_sub_sampled,Fitted_CTC,'b-');
    hold on;
    plot(time_vec_sub_sampled,Sim_Ct_low_res_noise,'k');
    hold off;
    
    
    if ( mod(i,10) == 0 )
        display('-I- Finished 10 iterations');
    end
    
    
end
