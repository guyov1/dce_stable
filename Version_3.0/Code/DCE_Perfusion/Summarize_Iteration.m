function [ results ] = Summarize_Iteration( Sim_Struct, Verbosity, iter_num, avg_num)

% Take from struct variables used in local function
Ktrans                                      = Sim_Struct.Ktrans;
Vb_larss                                = Sim_Struct.Vb_larss;
Ve_larss                                = Sim_Struct.Ve_larss;
additional_AIF_delay_sec                = Sim_Struct.additional_AIF_delay_sec;
est_sigma_noise_vec                     = Sim_Struct.est_sigma_noise_vec;
est_t_d_noise_vec                       = Sim_Struct.est_t_d_noise_vec;
est_amp_noise_vec                       = Sim_Struct.est_amp_noise_vec;
est_F_noise_vec                         = Sim_Struct.est_F_noise_vec;
est_Delay_sec_noise_vec                 = Sim_Struct.est_Delay_sec_noise_vec;
est_Delay_sec_using_Gaussian_noise_vec  = Sim_Struct.est_Delay_sec_using_Gaussian_noise_vec;
est_Ktrans_Patlak_noise_vec                 = Sim_Struct.est_Ktrans_Patlak_noise_vec;
est_Ktrans_Two_Comp_noise_vec               = Sim_Struct.est_Ktrans_Two_Comp_noise_vec;
est_E_noise_vec                         = Sim_Struct.est_E_noise_vec;
est_PS_noise_vec                        = Sim_Struct.est_PS_noise_vec;
est_Vb_Patlak_noise_vec                 = Sim_Struct.est_Vb_Patlak_noise_vec;
est_Vb_Two_Comp_noise_vec               = Sim_Struct.est_Vb_Two_Comp_noise_vec;
est_Vd_noise_vec                        = Sim_Struct.est_Vd_noise_vec;
est_Vd_normal_tis_noise_vec             = Sim_Struct.est_Vd_normal_tis_noise_vec;
est_Ve_Two_Comp_noise_vec               = Sim_Struct.est_Ve_Two_Comp_noise_vec;
est_MTT_noise_vec                       = Sim_Struct.est_MTT_noise_vec;
est_MTT_normal_tis_noise_vec            = Sim_Struct.est_MTT_normal_tis_noise_vec;
sigma                                   = Sim_Struct.sigma;
SNR_ratio                               = Sim_Struct.SNR_ratio;
sec_interval                            = Sim_Struct.sec_interval;
t_d                                     = Sim_Struct.t_d;
amplitude                               = Sim_Struct.amplitude;
F                                       = Sim_Struct.F;
E                                       = Sim_Struct.E;
PS                                      = Sim_Struct.PS;
Vd                                      = Sim_Struct.Vd;
MTT                                     = Sim_Struct.MTT;
num_averages                            = Sim_Struct.num_averages;
num_results_parameters                  = Sim_Struct.num_results_parameters;

est_F_Sourbron_noise_vec                = Sim_Struct.est_F_Sourbron_noise_vec;
est_Ktrans_Sourbron_Two_Comp_noise_vec      = Sim_Struct.est_Ktrans_Sourbron_Two_Comp_noise_vec;
est_Vb_Sourbron_Two_Comp_noise_vec      = Sim_Struct.est_Vb_Sourbron_Two_Comp_noise_vec;
est_Ve_Sourbron_Two_Comp_noise_vec      = Sim_Struct.est_Ve_Sourbron_Two_Comp_noise_vec;

if strcmp(Verbosity,'Full')
    display('-I- Summarizing iteration result...');
end

%% Calculate the mean value of all averages iterations
% Gaussian parameters
est_sigma_noise_avg                    = mean(est_sigma_noise_vec);
est_t_d_noise_avg                      = mean(est_t_d_noise_vec);
est_amp_noise_avg                      = mean(est_amp_noise_vec);
% Larsson Parameters
est_F_noise_avg                        = mean(est_F_noise_vec);
est_Delay_sec_noise_avg                = mean(est_Delay_sec_noise_vec);
est_Delay_sec_using_Gaussian_noise_avg = mean(est_Delay_sec_using_Gaussian_noise_vec);

est_Ktrans_Patlak_noise_avg      = mean(est_Ktrans_Patlak_noise_vec);
est_Ktrans_Two_Comp_noise_avg    = mean(est_Ktrans_Two_Comp_noise_vec);
est_E_noise_avg              = mean(est_E_noise_vec);
est_PS_noise_avg             = mean(est_PS_noise_vec);
est_Vb_Patlak_noise_avg      = mean(est_Vb_Patlak_noise_vec);
est_Vb_Two_Comp_noise_avg    = mean(est_Vb_Two_Comp_noise_vec);
est_Vd_noise_avg             = mean(est_Vd_noise_vec);
est_Vd_normal_tis_noise_avg  = mean(est_Vd_normal_tis_noise_vec);
est_MTT_noise_avg            = mean(est_MTT_noise_vec);
est_MTT_normal_tis_noise_avg = mean(est_MTT_normal_tis_noise_vec);
est_Ve_Two_Comp_noise_avg    = mean(est_Ve_Two_Comp_noise_vec);

est_F_Sourbron_noise_avg                = mean(est_F_Sourbron_noise_vec);
est_Ktrans_Sourbron_Two_Comp_noise_avg      = mean(est_Ktrans_Sourbron_Two_Comp_noise_vec);
est_Vb_Sourbron_Two_Comp_noise_avg      = mean(est_Vb_Sourbron_Two_Comp_noise_vec);
est_Ve_Sourbron_Two_Comp_noise_avg      = mean(est_Ve_Sourbron_Two_Comp_noise_vec);

% Calculate STD (standard deviations) for results
% Gaussian
est_sigma_noise_sec_vec      = 60 * est_sigma_noise_vec;
std_est_sigma                = std(est_sigma_noise_sec_vec);
std_est_t_d                  = std(est_t_d_noise_vec);
std_est_amp                  = std(est_amp_noise_vec);
% Larrson
std_est_F                    = std(est_F_noise_vec);
std_est_Delay                = std(est_Delay_sec_noise_vec);
std_est_Delay_using_Gaussian = std(est_Delay_sec_using_Gaussian_noise_vec);

std_Ktrans_Patlak                = std(est_Ktrans_Patlak_noise_vec);
std_Ktrans_Two_Comp              = std(est_Ktrans_Two_Comp_noise_vec);
std_E                        = std(est_E_noise_vec);
std_PS                       = std(est_PS_noise_vec);
std_Vb_Patlak                = std(est_Vb_Patlak_noise_vec);
std_Vb_Two_Comp              = std(est_Vb_Two_Comp_noise_vec);
std_Vd                       = std(est_Vd_noise_vec);
std_Vd_normal_tis            = std(est_Vd_normal_tis_noise_vec);
std_MTT                      = std(est_MTT_noise_vec);
std_MTT_normal_tis           = std(est_MTT_normal_tis_noise_vec);
std_Ve_Two_Comp              = std(est_Ve_Two_Comp_noise_vec);

% Sourbron
std_F_Sourbron               = std(est_F_Sourbron_noise_vec);
std_Ktrans_Sourbron              = std(est_Ktrans_Sourbron_Two_Comp_noise_vec);
std_Vb_Sourbron              = std(est_Vb_Sourbron_Two_Comp_noise_vec);
std_Ve_Sourbron              = std(est_Ve_Sourbron_Two_Comp_noise_vec);


% Collect results ( SNR, original params, estimated params)
% SNR, original sigma, estimated sigma, error percent, est. values STD
% original t_d, estimated t_d, error percent,  est. values STD

sigma_in_sec            = 60 * sigma(iter_num);
est_sigma_noise_avg_sec = 60 * est_sigma_noise_avg;

% Initiate results
results                 = zeros(num_results_parameters,1);

results(1)     = SNR_ratio(iter_num);
results(2)     = sec_interval(iter_num);

% Gaussian parameters
results(3)     = sigma_in_sec;
results(4)     = est_sigma_noise_avg_sec;
results(5)     = ( abs(sigma_in_sec - est_sigma_noise_avg_sec) / sigma_in_sec ) * 100;
results(6)     = std_est_sigma;

results(7)     = t_d(iter_num);
results(8)     = est_t_d_noise_avg;
results(9)     = ( abs(t_d(iter_num) - est_t_d_noise_avg) / t_d(iter_num) ) * 100;
results(10)    = std_est_t_d;

results(11)    = amplitude(iter_num);
results(12)    = est_amp_noise_avg;
results(13)    = ( abs(amplitude(iter_num) - est_amp_noise_avg) / amplitude(iter_num) ) * 100;
results(14)    = std_est_amp;

% Larsson parameters
results(15)    = F(iter_num);
results(16)    = est_F_noise_avg;
results(17)    = ( abs(F(iter_num) - est_F_noise_avg) / F(iter_num) ) * 100;
results(18)    = std_est_F;

results(19)    = additional_AIF_delay_sec(iter_num);
results(20)    = est_Delay_sec_noise_avg;
results(21)    = ( abs(additional_AIF_delay_sec(iter_num) - est_Delay_sec_noise_avg) / additional_AIF_delay_sec(iter_num) ) * 100;
results(22)    = std_est_Delay;

results(23)    = Ktrans(iter_num);
results(24)    = est_Ktrans_Patlak_noise_avg;
results(25)    = ( abs(Ktrans(iter_num) - est_Ktrans_Patlak_noise_avg) / Ktrans(iter_num) ) * 100;
results(26)    = std_Ktrans_Patlak;

results(27)    = Ktrans(iter_num);
results(28)    = est_Ktrans_Two_Comp_noise_avg;
results(29)    = ( abs(Ktrans(iter_num) - est_Ktrans_Two_Comp_noise_avg) / Ktrans(iter_num) ) * 100;
results(30)    = std_Ktrans_Two_Comp;

results(31)    = PS(iter_num);
results(32)    = est_PS_noise_avg;
results(33)    = ( abs(PS(iter_num) - est_PS_noise_avg) / PS(iter_num) ) * 100;
results(34)    = std_PS;

results(35)    = Vb_larss(iter_num);
results(36)    = est_Vb_Patlak_noise_avg;
results(37)    = ( abs(Vb_larss(iter_num) - est_Vb_Patlak_noise_avg) / Vb_larss(iter_num) ) * 100;
results(38)    = std_Vb_Patlak;

results(39)    = Vb_larss(iter_num);
results(40)    = est_Vb_Two_Comp_noise_avg;
results(41)    = ( abs(Vb_larss(iter_num) - est_Vb_Two_Comp_noise_avg) / Vb_larss(iter_num) ) * 100;
results(42)    = std_Vb_Two_Comp;

results(43)    = Vd(iter_num);
results(44)    = est_Vd_noise_avg;
results(45)    = ( abs(Vd(iter_num) - est_Vd_noise_avg) / Vd(iter_num) ) * 100;
results(46)    = std_Vd;

results(47)    = Vd(iter_num);
results(48)    = est_Vd_normal_tis_noise_avg;
results(49)    = ( abs(Vd(iter_num) - est_Vd_normal_tis_noise_avg) / Vd(iter_num) ) * 100;
results(50)    = std_Vd_normal_tis;

results(51)    = MTT(iter_num);
results(52)    = est_MTT_noise_avg;
results(53)    = ( abs(MTT(iter_num) - est_MTT_noise_avg) / MTT(iter_num) ) * 100;
results(54)    = std_MTT;

results(55)    = MTT(iter_num);
results(56)    = est_MTT_normal_tis_noise_avg;
results(57)    = ( abs(MTT(iter_num) - est_MTT_normal_tis_noise_avg) / MTT(iter_num) ) * 100;
results(58)    = std_MTT_normal_tis;

results(59)    = E(iter_num);
results(60)    = est_E_noise_avg;
results(61)    = ( abs(E(iter_num) - est_E_noise_avg) / E(iter_num) ) * 100;
results(62)    = std_E;

% AIF delay estimated using a gussian
results(63)    = additional_AIF_delay_sec(iter_num);
results(64)    = est_Delay_sec_using_Gaussian_noise_avg;
results(65)    = ( abs(additional_AIF_delay_sec(iter_num) - est_Delay_sec_using_Gaussian_noise_avg) / additional_AIF_delay_sec(iter_num) ) * 100;
results(66)    = std_est_Delay_using_Gaussian;

% Sourbron
results(67)    = F(iter_num);
results(68)    = est_F_Sourbron_noise_avg;
results(69)    = ( abs(F(iter_num) - est_F_Sourbron_noise_avg) / F(iter_num) ) * 100;
results(70)    = std_F_Sourbron;

results(71)    = Ktrans(iter_num);
results(72)    = est_Ktrans_Sourbron_Two_Comp_noise_avg;
results(73)    = ( abs(Ktrans(iter_num) - est_Ktrans_Sourbron_Two_Comp_noise_avg) / Ktrans(iter_num) ) * 100;
results(74)    = std_Ktrans_Sourbron;

results(75)    = Vb_larss(iter_num);
results(76)    = est_Vb_Sourbron_Two_Comp_noise_avg;
results(77)    = ( abs(Vb_larss(iter_num) - est_Vb_Sourbron_Two_Comp_noise_avg) / Vb_larss(iter_num) ) * 100;
results(78)    = std_Vb_Sourbron;

results(79)    = Ve_larss(iter_num);
results(80)    = est_Ve_Sourbron_Two_Comp_noise_avg;
results(81)    = ( abs(Ve_larss(iter_num) - est_Ve_Sourbron_Two_Comp_noise_avg) / Ve_larss(iter_num) ) * 100;
results(82)    = std_Ve_Sourbron;

% Continued Larsson paramters (Ve + absolute error

results(83)    = Ve_larss(iter_num);
results(84)    = est_Ve_Two_Comp_noise_avg;
results(85)    = ( abs(Ve_larss(iter_num) - est_Ve_Two_Comp_noise_avg) / Ve_larss(iter_num) ) * 100;
results(86)    = std_Ve_Two_Comp;


results(87)    = abs(F(iter_num)                        - est_F_noise_avg              );
results(88)    = abs(additional_AIF_delay_sec(iter_num) - est_Delay_sec_noise_avg      );
results(89)    = abs(Ktrans(iter_num)                       - est_Ktrans_Patlak_noise_avg      );
results(90)    = abs(Ktrans(iter_num)                       - est_Ktrans_Two_Comp_noise_avg    );
results(91)    = abs(PS(iter_num)                       - est_PS_noise_avg             );
results(92)    = abs(Vb_larss(iter_num)                 - est_Vb_Patlak_noise_avg      );
results(93)    = abs(Vb_larss(iter_num)                 - est_Vb_Two_Comp_noise_avg    );
results(94)    = abs(Vd(iter_num)                       - est_Vd_noise_avg             );
results(95)    = abs(Vd(iter_num)                       - est_Vd_normal_tis_noise_avg  );
results(96)    = abs(MTT(iter_num)                      - est_MTT_noise_avg            );
results(97)    = abs(MTT(iter_num)                      - est_MTT_normal_tis_noise_avg );
results(98)    = abs(E(iter_num)                        - est_E_noise_avg              );
results(99)    = abs(Ve_larss(iter_num)                 - est_Ve_Two_Comp_noise_avg    );

% Calculate iteration time
iter_finish = toc;

if ~strcmp(Verbosity,'None')
    display(sprintf('Finished iteration %d . Num Averages: %d. Took %.2f seconds to finish...',iter_num,num_averages,iter_finish));
end

end