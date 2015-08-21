function Filter = Adjusted_Larsson_Filter_no_Ve_no_E(time_vec_min, Fp, Vp)

% No extra-vascular - Uptake Model
Tp     = (Vp/Fp);
IRF    = Fp * exp(-(1/Tp) * time_vec_min);
Filter = IRF;

end