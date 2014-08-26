
function Filter = Adjusted_Larsson_Filter(time_vec_min, F, Vb, E, Ve)

chosen_time_stamps = length(time_vec_min);
chosen_time_vec    = time_vec_min;
Num_iterations     = 1;

% Create Larsson's filter (as described in article)
Vtis                   = ( 1 - (Vb./100) )*100; % Convert to [mL/100g]

PS     = ( E .* F ) ./ (1 - E);                             % [mL/100g/min]
alpha  = ( F + PS ) ./ Vb;                   % [1/min]
beta   = ( Vtis .* PS ) ./ (Vb.*Ve); % [1/min]
gamma  =  PS ./ Vtis;                         % [1/min]
theta  =  PS ./ Ve;               % [1/min]
a      = (1/2) * (theta + alpha + sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]
b      = (1/2) * (theta + alpha - sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]

IRF    = ( (a - theta - (PS/Vb) )*exp(-a*time_vec_min) - (b - theta - (PS/Vb) )*exp(-b*time_vec_min) ) / (a - b); % No units

% Vector implementation
%IRF    = ( repmat((a - theta - (PS./Vb) ),1,chosen_time_stamps) .* exp(-repmat(a,1,chosen_time_stamps).*repmat(chosen_time_vec,Num_iterations,1)) - ...
%           repmat((b - theta - (PS./Vb) ),1,chosen_time_stamps) .* exp(-repmat(b,1,chosen_time_stamps).*repmat(chosen_time_vec,Num_iterations,1)) ) ./ (repmat(a - b,1,chosen_time_stamps)); % No units

Filter = IRF;

end