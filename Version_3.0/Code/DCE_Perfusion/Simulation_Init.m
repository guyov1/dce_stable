function [ Sim_Struct ] = Simulation_Init( Sim_Struct, Verbosity )

% Take from struct variables used in local function
SNR_single                   = Sim_Struct.SNR_single;
num_iterations               = Sim_Struct.num_iterations;
additional_AIF_delay_sec_vec = Sim_Struct.additional_AIF_delay_sec_vec;
num_averages                 = Sim_Struct.num_averages;

if ~strcmp(Verbosity,'None')
    display('-I- Initiating iteration parameters...');
end

% Update SNR according to iteration number
if (Sim_Struct.iterate_SNR)
    % Each iteration decrease SNR
    Sim_Struct.SNR_ratio = Sim_Struct.SNR_vec;
else % Same SNR in each iteration
    Sim_Struct.SNR_ratio = repmat(SNR_single,1,num_iterations);
end

% Update time resolution according to iteration number
if (Sim_Struct.iterate_sec_interval)
    Sim_Struct.sec_interval   = Sim_Struct.sec_vec;     %[sec]
else
    Sim_Struct.sec_interval   = repmat(Sim_Struct.sec_interval, 1, num_iterations);     %[sec]
end

Sim_Struct.Fs             = 1 ./ Sim_Struct.sec_interval; %[Hz] - Sampling rate
Sim_Struct.min_interval   = Sim_Struct.sec_interval ./ 60;  %[min]

% Update AIF delay according to iteration number
if (Sim_Struct.iterate_AIF_delay)
    Sim_Struct.additional_AIF_delay_sec = additional_AIF_delay_sec_vec; % Delay added to AIF before filtering
    
else
    Sim_Struct.additional_AIF_delay_sec   = repmat(Sim_Struct.additional_AIF_delay_sec, 1, num_iterations);     %[sec]
end
Sim_Struct.additional_AIF_delay_min = Sim_Struct.additional_AIF_delay_sec/60; % Translate to minutes

% Change gaussian parameters in each iteration
if (Sim_Struct.iterate_gaussian_sigma)
    % Each iteration increase sigma
    %sigma = (2/60)  * iter_num;
    Sim_Struct.sigma = Sim_Struct.sigma_vec; % [min]
else
    Sim_Struct.sigma = repmat(Sim_Struct.sigma,1,num_iterations);
end
if (Sim_Struct.iterate_gaussian_time_delay)
    % Each iteration increase delay time
    Sim_Struct.t_d   = Sim_Struct.t_d_vec;  % Time shift in minutes
else
    Sim_Struct.t_d   = repmat(Sim_Struct.t_d,1,num_iterations);  % Time shift in minutes
end
if (Sim_Struct.iterate_gaussian_amplitude)
    % Each iteration increase amplitude
    Sim_Struct.amplitude = Sim_Struct.amplitude_vec;
else
    Sim_Struct.amplitude = repmat(Sim_Struct.amplitude,1,num_iterations);
end

% Create the the spline basis matrix

% Regular B-splines
Sim_Struct.B_mat      = Create_B_matrix(Sim_Struct.knots,Sim_Struct.time_vec_minutes,Sim_Struct.poly_deg-1);
% PCA Splines
Basis_Mat             = PCA_basis(Sim_Struct, Sim_Struct.time_vec_minutes);
% Take # of eigen-vectors similar to B-splines
num_cols_B_mat        = size(Sim_Struct.B_mat,2);
Sim_Struct.PCA_B_mat  = Basis_Mat(:,1:num_cols_B_mat);
% Laguerre splines
alpha = 0;
Sim_Struct.B_laguerre = Create_Laguerre_matrix(Sim_Struct.knots,Sim_Struct.time_vec_minutes,alpha,Sim_Struct.poly_deg-1);

% Plot PCA vs. Bsplines
% figure;
% subplot(1,2,1);
% plot(Sim_Struct.B_mat(:,1:5));
% xlabel('Sample #');
% ylabel('Amplitude');
% title('B-splines first 5 vectors');
% subplot(1,2,2);
% plot(Sim_Struct.PCA_B_mat(:,1:5));
% xlabel('Sample #');
% ylabel('Amplitude');
% title('PCA first 5 eigen-vectors');


% Gaussian parameters
Sim_Struct.est_sigma_noise_vec                    = zeros(1,num_averages);
Sim_Struct.est_t_d_noise_vec                      = zeros(1,num_averages);
Sim_Struct.est_amp_noise_vec                      = zeros(1,num_averages);

% Larsson's parameters
Sim_Struct.est_F_noise_vec                        = zeros(1,num_averages);
Sim_Struct.est_Delay_sec_noise_vec                = zeros(1,num_averages);
Sim_Struct.est_Delay_sec_using_Gaussian_noise_vec = zeros(1,num_averages);

Sim_Struct.est_Ktrans_Patlak_noise_vec                = zeros(1,num_averages);
Sim_Struct.est_Ktrans_Two_Comp_noise_vec              = zeros(1,num_averages);
Sim_Struct.est_E_noise_vec                        = zeros(1,num_averages);
Sim_Struct.est_PS_noise_vec                       = zeros(1,num_averages);
Sim_Struct.est_Vb_Patlak_noise_vec                = zeros(1,num_averages);
Sim_Struct.est_Vb_Two_Comp_noise_vec              = zeros(1,num_averages);
Sim_Struct.est_Ve_Two_Comp_noise_vec              = zeros(1,num_averages);
Sim_Struct.est_Vd_noise_vec                       = zeros(1,num_averages);
Sim_Struct.est_Vd_normal_tis_noise_vec            = zeros(1,num_averages);
Sim_Struct.est_MTT_noise_vec                      = zeros(1,num_averages);
Sim_Struct.est_MTT_normal_tis_noise_vec           = zeros(1,num_averages);

% Sourbron's parameters
Sim_Struct.est_F_Sourbron_noise_vec               = zeros(1,num_averages);
Sim_Struct.est_Ktrans_Sourbron_Two_Comp_noise_vec     = zeros(1,num_averages);
Sim_Struct.est_Vb_Sourbron_Two_Comp_noise_vec     = zeros(1,num_averages);
Sim_Struct.est_Ve_Sourbron_Two_Comp_noise_vec     = zeros(1,num_averages);

if strcmp(Verbosity,'Full')
    display('-I- Finished initiating iteration parameters...');
end

end