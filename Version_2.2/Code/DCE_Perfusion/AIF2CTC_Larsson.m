function [ CTC ] = AIF2CTC_Larsson( AIF, Flow, Vb, Ve, E, BAT, time_res_sec )
%AIF2CTC_Larsson
%   Creates CTC out of AIF and Larsson parameters
%   Input -
%           AIF          - Arterial input function
%           Flow         -
%           Vb           - intra-vascular volume
%           Ve           - extra-vascular extra-cellular volume
%           E            - Extraction fraction (permeability)
%           BAT          - Bolus arrival time of AIF
%           time_res_sec - AIF's time interval in seconds
%   Output -
%           CTC          - Concentration time curve


% Convert time to minutes
time_res_min    = time_res_sec/60;
% Number of time stamps
num_time_stamps = length(AIF);
chosen_time_vec = (1:num_time_stamps)*time_res_min;


x_dim      = size(Flow,1);
y_dim      = size(Flow,2);
z_dim      = size(Flow,3);

% Change 3D maps to 1-dim vector
F_vec      = reshape(Flow,x_dim*y_dim*z_dim,1);
Vb_vec     = reshape(Vb,x_dim*y_dim*z_dim,1);
Ve_vec     = reshape(Ve,x_dim*y_dim*z_dim,1);
E_vec      = reshape(E,x_dim*y_dim*z_dim,1);

num_filters = length(F_vec);

%% Create the filter out of parameters
Hct    = 0.38;
Ki     = E_vec .* F_vec;                             % [mL/100g/min]
alpha  = ( F_vec + Ki ) ./ Vb_vec;                   % [1/min]

Vtis   = ( 1 - (Vb_vec./100) )*100; % Convert to [mL/100g]

beta   = ( Vtis .* (1-Hct) .* Ki ) ./ (Vb_vec.*Ve_vec); % [1/min]
gamma  = Ki ./ Vtis;                         % [1/min]
theta  = (Ki*(1-Hct))   ./ Ve_vec;               % [1/min]
a      = (1/2) * (theta + alpha + sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]
b      = (1/2) * (theta + alpha - sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]

IRF    = ( repmat((a - theta - (Ki./Vb_vec) ),1,num_time_stamps) .* exp(-repmat(a,1,num_time_stamps).*repmat(chosen_time_vec,num_filters,1)) - ...
           repmat((b - theta - (Ki./Vb_vec) ),1,num_time_stamps) .* exp(-repmat(b,1,num_time_stamps).*repmat(chosen_time_vec,num_filters,1)) ) ./ (repmat(a - b,1,num_time_stamps)); % No units

Filter_high_res_matrix = repmat(F_vec,1,num_time_stamps) .* IRF;


%% Delay the AIF by BAT

AIF_shifted                     = lagmatrix(AIF,BAT);
AIF_shifted(isnan(AIF_shifted)) = 0;
AIF_shifted_T                   = AIF_shifted';

%% Pass the AIF through filters
CTC         = zeros(num_time_stamps,num_filters);
for i = 1:num_filters
    CTC(:,i) = filter(Filter_high_res_matrix(i,:)*time_res_min,1,AIF_shifted_T(i,:));
end


end

