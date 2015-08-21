function Filter = Adjusted_Larsson_Filter_no_E_High_F(time_vec_min, Vb)

% Highly Perfused and no indicator exchange

Vb_delta    = zeros(size(time_vec_min));
Vb_delta(1) = Vb;
IRF         = Vb_delta;


Filter = IRF;

end