
function Filter = Larsson_Filter(time_vec_min, F, Vb, E, Ve, Hct)

chosen_time_stamps = length(time_vec_min);
chosen_time_vec    = time_vec_min;
Num_iterations     = 1;

% Create Larsson's filter (as described in article)
Vtis                   = ( 1 - (Vb./100) )*100; % Convert to [mL/100g]

%Ki    = (1-Hct)*Ktrans;
Ki     =  E .* F;                             % [mL/100g/min]
alpha  = ( F + Ki ) ./ Vb;                   % [1/min]
beta   = ( Vtis .* (1-Hct) .* Ki ) ./ (Vb.*Ve); % [1/min]
gamma  =  Ki ./ Vtis;                         % [1/min]
theta  = (Ki.*(1-Hct))   ./ Ve;               % [1/min]
a      = (1/2) * (theta + alpha + sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]
b      = (1/2) * (theta + alpha - sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]

IRF    = ( (a - theta - (Ki/Vb) )*exp(-a*time_vec_min) - (b - theta - (Ki/Vb) )*exp(-b*time_vec_min) ) / (a - b); % No units

% Vector implementation
%IRF    = ( repmat((a - theta - (Ki./Vb) ),1,chosen_time_stamps) .* exp(-repmat(a,1,chosen_time_stamps).*repmat(chosen_time_vec,Num_iterations,1)) - ...
%           repmat((b - theta - (Ki./Vb) ),1,chosen_time_stamps) .* exp(-repmat(b,1,chosen_time_stamps).*repmat(chosen_time_vec,Num_iterations,1)) ) ./ (repmat(a - b,1,chosen_time_stamps)); % No units


Filter = IRF;

end