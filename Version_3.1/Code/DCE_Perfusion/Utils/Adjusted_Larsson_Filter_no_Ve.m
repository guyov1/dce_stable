function Filter = Adjusted_Larsson_Filter_no_Ve(time_vec_min, Fp, Vp, E)

% No extra-vascular - Uptake Model
Tp     = (Vp/Fp) * (1-E);
Ktrans = E*Fp;
IRF    = Fp * exp(-(1/Tp) * time_vec_min) + Ktrans*(1-exp(-(1/Tp) * time_vec_min));
Filter = IRF;

end