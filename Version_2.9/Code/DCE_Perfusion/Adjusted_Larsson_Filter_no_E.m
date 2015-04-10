function Filter = Adjusted_Larsson_Filter_no_E(time_vec_min, F, Vb)

% No indicator exchange

IRF = F * exp(-(F/Vb) * time_vec_min);

Filter = IRF;

end