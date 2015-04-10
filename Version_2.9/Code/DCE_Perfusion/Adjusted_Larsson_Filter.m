
function Filter = Adjusted_Larsson_Filter(time_vec_min, F, Vb, E, Ve)

%chosen_time_stamps = length(time_vec_min);
%chosen_time_vec    = time_vec_min;
%Num_iterations     = 1;

% In the , the IRF reduces to:
Vtis                   = ( 1 - (Vb./100) )*100; % Convert to [mL/100g]

PS     = ( E .* F ) ./ (1 - E);                             % [mL/100g/min]
alpha  = ( F + PS ) ./ Vb;                   % [1/min]
beta   = ( Vtis .* PS ) ./ (Vb.*Ve); % [1/min]

gamma  =  PS ./ Vtis;                         % [1/min]
theta  =  PS ./ Ve;               % [1/min]
a      = (1/2) * (theta + alpha + sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]
b      = (1/2) * (theta + alpha - sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]

IRF    = ( (a - theta - (PS/Vb) )*exp(-a*time_vec_min) - (b - theta - (PS/Vb) )*exp(-b*time_vec_min) ) / (a - b); % No units



% if ( F > exp(20) ) % Highly Perfused
%
%     PS          = E; % In this case, whould make sure that E is Ktrans (and also PS = Ktrans)
%     Vb_delta    = zeros(size(time_vec_min));
%     Vb_delta(1) = Vb;
%     if (E == 0 || Ve ==0)
%         IRF         = Vb_delta;
%     else
%         IRF         = Vb_delta + ( PS * exp(-(PS/Ve) * time_vec_min) );
%     end
%
% elseif (E == 0)    % No indicator exchange model
%     IRF = F * exp(-(F/Vb) * time_vec_min);
% else
%     % Create Larsson's filter (as described in article)
%     Vtis                   = ( 1 - (Vb./100) )*100; % Convert to [mL/100g]
%
%     PS     = ( E .* F ) ./ (1 - E);                             % [mL/100g/min]
%     alpha  = ( F + PS ) ./ Vb;                   % [1/min]
%     beta   = ( Vtis .* PS ) ./ (Vb.*Ve); % [1/min]
%
%     gamma  =  PS ./ Vtis;                         % [1/min]
%     theta  =  PS ./ Ve;               % [1/min]
%     a      = (1/2) * (theta + alpha + sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]
%     b      = (1/2) * (theta + alpha - sqrt(theta.^2 + alpha.^2 - 2*theta.*alpha + 4*gamma.*beta) ); % [1/min]
%
%     IRF    = ( (a - theta - (PS/Vb) )*exp(-a*time_vec_min) - (b - theta - (PS/Vb) )*exp(-b*time_vec_min) ) / (a - b); % No units
%
%     % Vector implementation
%     %IRF    = ( repmat((a - theta - (PS./Vb) ),1,chosen_time_stamps) .* exp(-repmat(a,1,chosen_time_stamps).*repmat(chosen_time_vec,Num_iterations,1)) - ...
%     %           repmat((b - theta - (PS./Vb) ),1,chosen_time_stamps) .* exp(-repmat(b,1,chosen_time_stamps).*repmat(chosen_time_vec,Num_iterations,1)) ) ./ (repmat(a - b,1,chosen_time_stamps)); % No units
% end



Filter = IRF;

end