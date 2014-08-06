% Uses PCA to get the best possible basis functions for the bi-exp filter

function [] = PCA_basis(sec_inter)

display(['-I- Starting PCA basis creation for:' sec_inter ' seconds interval...']);

%% Filter generation

% Larsson filter parameters
F_low_val  = 0;
F_max_val  = 140;
Vb_low_val = 3;
Vb_max_val = 20;
Ve_low_val = 1;
Ve_max_val = 10;
E_low_val  = 0;
E_max_val  = 1;

Hct        = 0.38;

Num_iterations = 10000;

F_vec      = F_low_val  + (F_max_val  - F_low_val) * rand(Num_iterations,1);
Vb_vec     = Vb_low_val + (Vb_max_val - Vb_low_val)* rand(Num_iterations,1);
Ve_vec     = Ve_low_val + (Ve_max_val - Vb_low_val)* rand(Num_iterations,1);
E_vec      = E_low_val  + (E_max_val  - E_low_val) * rand(Num_iterations,1);

% Time parameters
sec_interval_low_res     = sec_inter;
sec_interval_high_res    = 0.01;
min_interval_low_res     = sec_interval_low_res/60;
min_interval_high_res    = sec_interval_high_res/60;
total_sim_time_min       = 4;
num_time_stamps_high_res = round(total_sim_time_min/min_interval_high_res);
num_time_stamps_low_res  = round(total_sim_time_min/min_interval_low_res);
time_vec_high_resolution = (0:num_time_stamps_high_res-1) * min_interval_high_res;
time_vec_sub_sampled     = (0:num_time_stamps_low_res-1) * min_interval_low_res;
UpSampFactor             = round(min_interval_low_res / min_interval_high_res) ;

% Create Larsson's filter (as described in article)
chosen_time_vec        = time_vec_sub_sampled;
chosen_time_stamps     = num_time_stamps_low_res;
Vtis                   = ( 1 - (Vb_vec./100) )*100; % Convert to [mL/100g]

%Ki    = (1-Hct)*Ktrans;
Ki     = E_vec .* F_vec;                             % [mL/100g/min]
alpha  = ( F_vec + Ki ) ./ Vb_vec;                   % [1/min]

beta   = ( Vtis .* (1-Hct) .* Ki ) ./ (Vb_vec.*Ve_vec); % [1/min]
gamma  = Ki ./ Vtis;                             % [1/min]
theta  = (Ki*(1-Hct))   ./ Ve_vec;               % [1/min]
a      = (1/2) * (theta + alpha + sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]
b      = (1/2) * (theta + alpha - sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]

IRF    = ( repmat((a - theta - (Ki./Vb_vec) ),1,chosen_time_stamps) .* exp(-repmat(a,1,chosen_time_stamps).*repmat(chosen_time_vec,Num_iterations,1)) - ...
    repmat((b - theta - (Ki./Vb_vec) ),1,chosen_time_stamps) .* exp(-repmat(b,1,chosen_time_stamps).*repmat(chosen_time_vec,Num_iterations,1)) ) ./ (repmat(a - b,1,chosen_time_stamps)); % No units

Filter_high_res_matrix = repmat(F_vec,1,chosen_time_stamps) .* IRF;

%% Apllying PCA

[COEFF,SCORE,latent] = princomp(Filter_high_res_matrix); % returns latent, a vector containing the eigenvalues of the covariance matrix of X.

% Get the first vectors which explains 99% of observations
Ratio_explained = cumsum(latent ./ sum(latent));
[num_eigen_vecs] = find(Ratio_explained>0.99,1);

figure;
plot(chosen_time_vec,COEFF(:,1:num_eigen_vecs));
title('Eigen Vectors - Over 99% Var');
B_check = COEFF;
file_name = ['\\fmri-guy2\Dropbox\University\Msc\Thesis\General\Matlab Simulations\Flow Extraction\Latest Code\Utils\Alt_Basis_' num2str(sec_interval_low_res) '_sec' '.mat'];
save(file_name,'B_check');

% % Normalize the examples
% [Filter_high_res_matrix_norm, mu, sigma] = featureNormalize(Filter_high_res_matrix);
%
% [U, S] = pca(Filter_high_res_matrix_norm);
%
% %  Project the data onto K = num_time_stamps_low_res dimension
% K = num_time_stamps_low_res;
% Z = projectData(X_norm, U, K);
%
end