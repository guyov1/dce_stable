% Uses PCA to get the best possible basis functions for the bi-exp filter
function [B_PCA] = PCA_basis(Sim_Struct, time_vec)

display('-I- Creating PCA basis matrix...');

%% Get relevant parameters from simulation struct
% Time parameters
% sec_inter             = Sim_Struct.sec_interval(1);
% total_sim_time        = Sim_Struct.total_sim_time_min;
% min_interval_high_res = Sim_Struct.High_res_min;
% Larsson filter parameters
F_low                 = Sim_Struct.F_low  ;
F_max                 = Sim_Struct.F_max  ;
Vb_low                = Sim_Struct.Vb_low ;
Vb_max                = Sim_Struct.Vb_max ;
Ve_low                = Sim_Struct.Ve_low ;
Ve_max                = Sim_Struct.Ve_max ;
E_low                 = Sim_Struct.E_low  ;
E_max                 = Sim_Struct.E_max  ;
% More parameters
Num_iterations        = Sim_Struct.Num_iterations_PCA;
Hct                   = Sim_Struct.Hct_single;
adjusted_larsson      = Sim_Struct.Adjusted_Larsson_Model;

%% Filter generation
F      = F_low  + (F_max  - F_low) * rand(Num_iterations,1);
Vb     = Vb_low + (Vb_max - Vb_low)* rand(Num_iterations,1);
Ve     = Ve_low + (Ve_max - Vb_low)* rand(Num_iterations,1);
E      = E_low  + (E_max  - E_low) * rand(Num_iterations,1);

% Time parameters
% sec_interval_low_res     = sec_inter;
% min_interval_low_res     = sec_interval_low_res/60;
% 
% num_time_stamps_high_res = round(total_sim_time/min_interval_high_res);
% time_vec_high_resolution = (0 : num_time_stamps_high_res - 1) * min_interval_high_res;
% UpSampFactor             = round(min_interval_low_res / min_interval_high_res) ;
% 
% num_time_stamps_low_res  = round(total_sim_time/min_interval_low_res);
% time_vec_sub_sampled     = (0 : num_time_stamps_low_res  - 1) * min_interval_low_res;
% % Create Larsson's filter (as described in article)
% time_vec          = time_vec_sub_sampled;
% num_time_stamps   = num_time_stamps_low_res;
 
num_time_stamps = length(time_vec);

if adjusted_larsson
    Vtis                   = ( 1 - (Vb./100) )*100; % Convert to [mL/100g]
    
    PS     = ( E .* F ) ./ (1 - E);                             % [mL/100g/min]
    alpha  = ( F + PS ) ./ Vb;                   % [1/min]
    beta   = ( Vtis .* PS ) ./ (Vb.*Ve); % [1/min]
    gamma  =  PS ./ Vtis;                         % [1/min]
    theta  =  PS ./ Ve;               % [1/min]
    a      = (1/2) * (theta + alpha + sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]
    b      = (1/2) * (theta + alpha - sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]
    
    %IRF    = ( (a - theta - (PS/Vb) )*exp(-a*time_vec_min) - (b - theta - (PS/Vb) )*exp(-b*time_vec_min) ) / (a - b); % No units
    
    IRF    = ( repmat((a - theta - (PS./Vb) ),1,num_time_stamps) .* exp(-repmat(a,1,num_time_stamps).*repmat(time_vec,Num_iterations,1)) - ...
               repmat((b - theta - (PS./Vb) ),1,num_time_stamps) .* exp(-repmat(b,1,num_time_stamps).*repmat(time_vec,Num_iterations,1)) ) ./ (repmat(a - b,1,num_time_stamps)); % No units
else
    % Create Larsson's filter (as described in article)
    Vtis   = ( 1 - (Vb./100) )*100; % Convert to [mL/100g]
    
    %Ktrans    = (1-Hct)*Ktrans;
    Ktrans     =  E .* F;                             % [mL/100g/min]
    alpha  = ( F + Ktrans ) ./ Vb;                   % [1/min]
    beta   = ( Vtis .* (1-Hct) .* Ktrans ) ./ (Vb.*Ve); % [1/min]
    gamma  =  Ktrans ./ Vtis;                         % [1/min]
    theta  = (Ktrans*(1-Hct))   ./ Ve;               % [1/min]
    a      = (1/2) * (theta + alpha + sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]
    b      = (1/2) * (theta + alpha - sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]
    
    %IRF    = ( (a - theta - (Ktrans/Vb) )*exp(-a*time_vec_min) - (b - theta - (Ktrans/Vb) )*exp(-b*time_vec_min) ) / (a - b); % No units
    
    IRF    = ( repmat((a - theta - (Ktrans./Vb) ),1,num_time_stamps) .* exp(-repmat(a,1,num_time_stamps).*repmat(time_vec,Num_iterations,1)) - ...
        repmat((b - theta - (Ktrans./Vb) ),1,num_time_stamps) .* exp(-repmat(b,1,num_time_stamps).*repmat(time_vec,Num_iterations,1)) ) ./ (repmat(a - b,1,num_time_stamps)); % No units
end

Filter_high_res_matrix = repmat(F,1,num_time_stamps) .* IRF;

%% Apllying PCA

[COEFF,SCORE,latent] = princomp(Filter_high_res_matrix); % returns latent, a vector containing the eigenvalues of the covariance matrix of X.

% Get the first vectors which explains 99% of observations
Ratio_explained  = cumsum(latent ./ sum(latent));
[num_eigen_vecs] = find(Ratio_explained > 0.99,1);

% figure;
% plot(time_vec,COEFF(:,1:num_eigen_vecs));
% title('Eigen Vectors - Over 99% Var');

B_PCA     = COEFF;

%file_name = ['\\fmri-guy2\Dropbox\University\Msc\Thesis\General\Matlab Simulations\Flow Extraction\Latest Code\Utils\Alt_Basis_' num2str(sec_interval_low_res) '_sec_Total_' num2str(Sim_Struct.total_sim_time_min) '_min' '.mat'];
%save(file_name,'B_PCA');

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