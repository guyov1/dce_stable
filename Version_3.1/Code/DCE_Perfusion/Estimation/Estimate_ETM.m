function [ Sim_Struct ] = Estimate_ETM( Sim_Struct, Verbosity )

if ~strcmp(Verbosity,'None')
    display('-I- Starting Murase/Tofts Comparison...');
end

% Take from struct variables used in local function
min_interval                       = Sim_Struct.min_interval;
time_vec_minutes                   = Sim_Struct.time_vec_minutes;
Sim_AIF_delayed_no_noise           = Sim_Struct.Sim_AIF_delayed_no_noise;
Sim_AIF_with_noise                 = Sim_Struct.Sim_AIF_with_noise;
SNR_ratio                          = Sim_Struct.SNR_ratio;
One_Iteration_Murase_Tofts         = Sim_Struct.One_Iteration_Murase_Tofts;
Iterate_Murase_Tofts_num_iter      = Sim_Struct.Iterate_Murase_Tofts_num_iter;
algorithm_options                  = Sim_Struct.algorithm_options;
adjusted_larsson                   = Sim_Struct.Adjusted_Larsson_Model;
num_iterations                     = Sim_Struct.num_iterations;
Sim_Ct_larss_Murase_noise          = Sim_Struct.Sim_Ct_larss_kernel_noise;
Sim_Ct_larss_Murase_noise_high_res = Sim_Struct.Sim_Ct_larss_kernel_noise_high_res;

% Check if only one iteration is wanted
if (One_Iteration_Murase_Tofts)
    Iterate_Murase_Tofts_num_iter = 1;
end

% Initiate parameter matrices
Murase_params    = zeros(3,num_iterations);

% Initiate time calculation for first index
tic;
display(sprintf('-I- Starting simulation for %d voxels...',Iterate_Murase_Tofts_num_iter));

%parfor idx = 1 : Iterate_Murase_Tofts_num_iter
for idx = 1 : num_iterations
    
    % Create vectors of A matrix
    A_1 =  cumtrapz(time_vec_minutes,Sim_AIF_with_noise(:,idx)');
    A_2 = -cumtrapz(time_vec_minutes,Sim_Ct_larss_Murase_noise(:,idx)');
    A_3 =  Sim_AIF_with_noise(:,idx)';
    A   =  [A_1' A_2' A_3'];
   
    % Create c vector
    C_vec = Sim_Ct_larss_Murase_noise(:,idx);
    
    B     = A \ C_vec;
    %B = pinv(A) * C_vec;
    
    K_trans_Tofts_Murase = B(1) - B(2)*B(3);
    K_ep_Tofts_Murase    = B(2);
    Vp_Tofts_Murase      = B(3);
    
    % Make sure result is not negative
    K_trans_Tofts_Murase = max(0,K_trans_Tofts_Murase);
    K_ep_Tofts_Murase    = max(0,K_ep_Tofts_Murase);
    Vp_Tofts_Murase      = max(0,Vp_Tofts_Murase);
    
    Murase_params(:,idx)    = [K_trans_Tofts_Murase K_ep_Tofts_Murase Vp_Tofts_Murase];

end

display(sprintf('Finished simulation for %d voxels...',Iterate_Murase_Tofts_num_iter));
time_finish = toc;
display(sprintf('Took %.2f seconds to finish...',time_finish));
tic;

Sim_Struct.Est_Ktrans_vec = Murase_params(1,:);
Sim_Struct.Est_Kep_vec    = Murase_params(2,:);
Sim_Struct.Est_Vp_vec     = Murase_params(3,:);
Sim_Struct.Est_Ve_vec     = Sim_Struct.Est_Ktrans_vec ./ Sim_Struct.Est_Kep_vec;

%% Plot
% Plot a single iteration
plot_idx                   = Sim_Struct.ETM_idx_to_plot;
fixed_vp_vec_just_const    = [B(1)-B(2)*B(3) ; B(2); Sim_Struct.Vp_ETM(plot_idx)];
fixed_vp_vec               = [B(1)-B(2)*Sim_Struct.Vp_ETM(plot_idx) ; B(2); Sim_Struct.Vp_ETM(plot_idx)];

true_params_vec            = [Sim_Struct.Vp_ETM(plot_idx) Sim_Struct.Ktrans_ETM(plot_idx) Sim_Struct.kep_ETM(plot_idx)];
est_params_vec             = [Murase_params(3,plot_idx) Murase_params(1,plot_idx) Murase_params(2,plot_idx)];

figure;
subplot(2,1,1);
h1 = plot(Sim_Struct.time_vec_minutes_high_res,Sim_Ct_larss_Murase_noise_high_res(:,plot_idx)','--k');
hold on;
h2 = plot(time_vec_minutes,Sim_Ct_larss_Murase_noise(:,plot_idx)','*g');
h3 = plot(time_vec_minutes,A*B,'or');
h4 = plot(time_vec_minutes,A*fixed_vp_vec,'xb');
h5 = plot(time_vec_minutes,A*fixed_vp_vec_just_const,'dc');
hold off;

legend([h1 h2 h3 h4 h5], 'Cont. CTC','Sampled CTC', 'Estimated Murase Fit', 'Murase - Fit with adjusted Vp to real data', 'Murase - Fit with adjusted Vp to real data - Just Const');
title(['Murase Fit. True Params(Vp,Kt,Kep):  ' num2str(true_params_vec) ' . Estimed:  ' num2str(est_params_vec) ]);

subplot(2,1,2);


filter_true = ETM_Filter(Sim_Struct.time_vec_minutes_high_res, Sim_Struct.Vp_ETM(plot_idx), Sim_Struct.Ktrans_ETM(plot_idx), Sim_Struct.kep_ETM(plot_idx));
filter_est  = ETM_Filter(Sim_Struct.time_vec_minutes_high_res, Sim_Struct.Est_Vp_vec(plot_idx), Sim_Struct.Est_Ktrans_vec(plot_idx), Sim_Struct.Est_Kep_vec(plot_idx));
filter_est_fixed_vp  = ETM_Filter(Sim_Struct.time_vec_minutes_high_res, fixed_vp_vec(3), fixed_vp_vec(1), fixed_vp_vec(2));
filter_est_fixed_vp_just_const  = ETM_Filter(Sim_Struct.time_vec_minutes_high_res, fixed_vp_vec_just_const(3), fixed_vp_vec_just_const(1), fixed_vp_vec_just_const(2));

relevant_AIF = Sim_Struct.Sim_AIF_HighRes_delayed_no_noise(:,plot_idx);
CTC_true     = filter(filter_true*Sim_Struct.High_res_min,1,relevant_AIF) + Sim_Struct.Vp_ETM(plot_idx)*relevant_AIF;
CTC_est      = filter(filter_est*Sim_Struct.High_res_min,1,relevant_AIF) + Sim_Struct.Est_Vp_vec(plot_idx)*relevant_AIF;
CTC_est_fixed_vp            = filter(filter_est_fixed_vp*Sim_Struct.High_res_min,1,relevant_AIF) + fixed_vp_vec(3)*relevant_AIF;
CTC_est_fixed_vp_just_const = filter(filter_est_fixed_vp_just_const*Sim_Struct.High_res_min,1,relevant_AIF) + fixed_vp_vec_just_const(3)*relevant_AIF;


h6 = plot(Sim_Struct.time_vec_minutes_high_res,CTC_est,'*g');
hold on;
h7 = plot(Sim_Struct.time_vec_minutes_high_res,CTC_true,'or');
h8 = plot(Sim_Struct.time_vec_minutes_high_res,CTC_est_fixed_vp,'ms');
h9 = plot(Sim_Struct.time_vec_minutes_high_res,CTC_est_fixed_vp_just_const,'*y');
hold off;
legend([h6 h7 h8 h9],'True CTC','Est. CTC','Est. CTC fixed Vp','Est. CTC fixed Vp - Just Const'); 
title('Filtering process');







if strcmp(Verbosity,'Full')
    display('-I- Finished Murase/Tofts Comparison...');
end

end