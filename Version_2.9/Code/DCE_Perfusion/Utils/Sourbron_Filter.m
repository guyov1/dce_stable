function [ IRF_sourbron ] = Sourbron_Filter(time_vec_minutes, Fp, Vp, Fe, Ve)
% Create Sourbron's filter (as described in article)

T_P = Vp / (Fe + Fp);
T_E = Ve / Fe;
T_B = Vp / Fp; 

K_plus  = (1/2) * ( (T_P)^(-1) + (T_E)^(-1) + sqrt( ( (T_P)^(-1) + (T_E)^(-1) )^2 - ( 4 /(T_E*T_B) ) ) );
K_minus = (1/2) * ( (T_P)^(-1) + (T_E)^(-1) - sqrt( ( (T_P)^(-1) + (T_E)^(-1) )^2 - ( 4 /(T_E*T_B) ) ) );

E_minus      = ( K_plus - (T_B)^(-1) ) / ( K_plus - K_minus );

IRF_sourbron = exp(-time_vec_minutes*K_plus) + E_minus*(exp(-time_vec_minutes*K_minus)-exp(-time_vec_minutes*K_plus));


end

