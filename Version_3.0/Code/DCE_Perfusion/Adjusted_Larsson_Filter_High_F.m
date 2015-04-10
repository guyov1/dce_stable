
function Filter = Adjusted_Larsson_Filter_High_F(time_vec_min, Vb, Ktrans, Ve)

% Highly Perfused

% In this case, PS = Ktrans

Vb_delta    = zeros(size(time_vec_min));
Vb_delta(1) = Vb;
IRF         = Vb_delta + ( Ktrans * exp(-(Ktrans/Ve) * time_vec_min) );



Filter = IRF;

end