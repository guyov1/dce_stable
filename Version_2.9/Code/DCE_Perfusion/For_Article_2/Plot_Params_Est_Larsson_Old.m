


% load(DataPath_1);
% 
% 
% plot_error_results_flag_1           = plot_error_results_flag;
% iterate_SNR_1                       = iterate_SNR;                
% iterate_sec_interval_1              = iterate_sec_interval;       
% Ignore_Gaussian_Calculation_1       = Ignore_Gaussian_Calculation;
% iterate_gaussian_sigma_1            = iterate_gaussian_sigma;     
% iterate_gaussian_time_delay_1       = iterate_gaussian_time_delay;
% iterate_gaussian_amplitude_1        = iterate_gaussian_amplitude; 
% iterate_F_larsson_1                 = iterate_F_larsson;          
% iterate_Vb_larsson_1                = iterate_Vb_larsson;         
% iterate_E_larsson_1                 = iterate_E_larsson;           
% Vb_single_1                         = Vb_single;                   
% E_single_1                          = E_single;                  
% F_single_1                          = F_single;                      
% iterate_AIF_delay_1                 = iterate_AIF_delay;              
% iterate_uniformly_1                 = iterate_uniformly;             
% results_1                           = results;                      
% Filter_Est_Chosen_1                 = Filter_Est_Chosen;                 
% F_low_1                             = F_low;                        
% E_low_1                             = E_low;                        
% Vb_low_1                            = Vb_low;                         
% Ve_low_1                            = Ve_low;                       
% F_max_1                             = F_max;                      
% E_max_1                             = E_max;                         
% Vb_max_1                            = Vb_max;                         
% Ve_max_1                            = Ve_max;                         
% 
% % F
% real_larsson_F_vec_1          = results_1(15,:);
% est_larsson_F_vec_1           = results_1(16,:);
% error_percent_F_1             = results_1(17,:);
% std_F_1                       = results_1(18,:);
% 
% % AIF Delay time (Larsson)
% real_t_d_Larss_vec_sec_1      = results_1(19,:);
% est_t_d_Larss_vec_sec_1       = results_1(20,:);
% error_percent_t_d_Larss_1     = results_1(21,:);
% std_t_d_Larss_sec_1           = results_1(22,:);
% 
% real_larsson_Ktrans_2CXM_vec_1    = results_1(27,:);
% est_larsson_Ktrans_2CXM_vec_1     = results_1(28,:);
% error_percent_Ktrans_2CXM_1       = results_1(29,:);
% std_Ktrans_2CXM_1                 = results_1(30,:);
% 
% real_larsson_Vb_Patlak_vec_1  = results_1(35,:);
% est_larsson_Vb_Patlak_vec_1   = results_1(36,:);
% error_percent_Vb_Patlak_1     = results_1(37,:);
% std_Vb_Patlak_1               = results_1(38,:);
% real_larsson_Vb_2CXM_vec_1    = results_1(39,:);
% est_larsson_Vb_2CXM_vec_1     = results_1(40,:);
% error_percent_Vb_2CXM_1       = results_1(41,:);
% std_Vb_2CXM_1                 = results_1(42,:);
% 
% % E
% real_larsson_E_vec_1          = results_1(59,:);
% est_larsson_E_vec_1           = results_1(60,:);
% error_percent_E_1             = results_1(61,:);
% std_E_1                       = results_1(62,:);
% 
% % Larsson - Ve
% real_larsson_Ve_2CXM_vec_1    = results_1(83,:);
% est_larsson_Ve_2CXM_vec_1     = results_1(84,:);
% error_percent_Ve_2CXM_1       = results_1(85,:);
% std_Ve_2CXM_1                 = results_1(86,:);
% 
% %% Load data 2
% 
% load(DataPath_2);
% 
% 
% plot_error_results_flag_2           = plot_error_results_flag;
% iterate_SNR_2                       = iterate_SNR;
% iterate_sec_interval_2              = iterate_sec_interval;
% Ignore_Gaussian_Calculation_2       = Ignore_Gaussian_Calculation;
% iterate_gaussian_sigma_2            = iterate_gaussian_sigma;
% iterate_gaussian_time_delay_2       = iterate_gaussian_time_delay;
% iterate_gaussian_amplitude_2        = iterate_gaussian_amplitude;
% iterate_F_larsson_2                 = iterate_F_larsson;
% iterate_Vb_larsson_2                = iterate_Vb_larsson;
% iterate_E_larsson_2                 = iterate_E_larsson;
% Vb_single_2                         = Vb_single;
% E_single_2                          = E_single;
% F_single_2                          = F_single;
% iterate_AIF_delay_2                 = iterate_AIF_delay;
% iterate_uniformly_2                 = iterate_uniformly;
% results_2                           = results;
% Filter_Est_Chosen_2                 = Filter_Est_Chosen;
% F_low_2                             = F_low;
% E_low_2                             = E_low;
% Vb_low_2                            = Vb_low;
% Ve_low_2                            = Ve_low;
% F_max_2                             = F_max;
% E_max_2                             = E_max;
% Vb_max_2                            = Vb_max;
% Ve_max_2                            = Ve_max;
% 
% % F
% real_larsson_F_vec_2          = results_2(15,:);
% est_larsson_F_vec_2           = results_2(16,:);
% error_percent_F_2             = results_2(17,:);
% std_F_2                       = results_2(18,:);
% 
% % AIF Delay time (Larsson)
% real_t_d_Larss_vec_sec_2      = results_2(19,:);
% est_t_d_Larss_vec_sec_2       = results_2(20,:);
% error_percent_t_d_Larss_2     = results_2(21,:);
% std_t_d_Larss_sec_2           = results_2(22,:);
% 
% real_larsson_Ktrans_2CXM_vec_2    = results_2(27,:);
% est_larsson_Ktrans_2CXM_vec_2     = results_2(28,:);
% error_percent_Ktrans_2CXM_2       = results_2(29,:);
% std_Ktrans_2CXM_2                 = results_2(30,:);
% 
% real_larsson_Vb_Patlak_vec_2  = results_2(35,:);
% est_larsson_Vb_Patlak_vec_2   = results_2(36,:);
% error_percent_Vb_Patlak_2     = results_2(37,:);
% std_Vb_Patlak_2               = results_2(38,:);
% real_larsson_Vb_2CXM_vec_2    = results_2(39,:);
% est_larsson_Vb_2CXM_vec_2     = results_2(40,:);
% error_percent_Vb_2CXM_2       = results_2(41,:);
% std_Vb_2CXM_2                 = results_2(42,:);
% 
% % E
% real_larsson_E_vec_2          = results_2(59,:);
% est_larsson_E_vec_2           = results_2(60,:);
% error_percent_E_2             = results_2(61,:);
% std_E_2                       = results_2(62,:);
% 
% % Larsson - Ve
% real_larsson_Ve_2CXM_vec_2    = results_2(83,:);
% est_larsson_Ve_2CXM_vec_2     = results_2(84,:);
% error_percent_Ve_2CXM_2       = results_2(85,:);
% std_Ve_2CXM_2                 = results_2(86,:);
