function Filter = ETM_Filter(time_vec_min, Vp, Ktrans, kep)

%chosen_time_stamps = length(time_vec_min);
%chosen_time_vec    = time_vec_min;
%Num_iterations     = 1;

% % In the , the IRF reduces to:
% Vtis                   = ( 1 - (Vb./100) )*100; % Convert to [mL/100g]
% 
% PS     = ( E .* F ) ./ (1 - E);                             % [mL/100g/min]
% alpha  = ( F + PS ) ./ Vb;                   % [1/min]
% beta   = ( Vtis .* PS ) ./ (Vb.*Ve); % [1/min]
% 
% gamma  =  PS ./ Vtis;                         % [1/min]
% theta  =  PS ./ Ve;               % [1/min]
% a      = (1/2) * (theta + alpha + sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]
% b      = (1/2) * (theta + alpha - sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]

%kep      = Ktrans ./ Ve;
%IRF    = Ktrans*exp(-kep*time_vec_min); % No units

% Vp_delta = zeros(size(time_vec_min));
% Vp_delta(1) = Vp;
% Vp_delta = Vp * dirac(time_vec_min);
% IRF      = Vp_delta + Ktrans*exp(-kep*time_vec_min); % No units

% Can't handle dirac, use just seconds part
IRF      = Ktrans*exp(-kep*time_vec_min); % No units

Filter = IRF;

end