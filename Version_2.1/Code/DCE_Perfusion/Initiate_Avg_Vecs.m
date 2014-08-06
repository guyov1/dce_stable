% Initiate local loop variabels

% Calculate iteration time
tic;

% Initiate local loop vectors
num_averages                             = Sim_Struct_Replicated(iter_num).num_averages;
est_sigma_noise_vec                      = zeros(num_averages,1);
est_t_d_noise_vec                        = zeros(num_averages,1);
est_amp_noise_vec                        = zeros(num_averages,1);
est_F_noise_vec                          = zeros(num_averages,1);
est_Delay_sec_noise_vec                  = zeros(num_averages,1);
est_Delay_sec_using_Gaussian_noise_vec   = zeros(num_averages,1);
est_Ki_Patlak_noise_vec                  = zeros(num_averages,1);
est_Ki_Two_Comp_noise_vec                = zeros(num_averages,1);
est_Vb_Patlak_noise_vec                  = zeros(num_averages,1);
est_E_noise_vec                          = zeros(num_averages,1);
est_PS_noise_vec                         = zeros(num_averages,1);
est_Vb_Two_Comp_noise_vec                = zeros(num_averages,1);
est_Vd_noise_vec                         = zeros(num_averages,1);
est_Vd_normal_tis_noise_vec              = zeros(num_averages,1);
est_MTT_noise_vec                        = zeros(num_averages,1);
est_MTT_normal_tis_noise_vec             = zeros(num_averages,1);