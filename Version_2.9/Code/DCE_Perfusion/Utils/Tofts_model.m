
function to_fit = Tofts_model(time_vec_min, min_interval, AIF, Ktrans, Kep ,Vp)
   
    AIF_part        = Vp*AIF;
    Kep_Filter_Part = filter(Ktrans*exp(-Kep*time_vec_min)*min_interval,1,AIF);
    to_fit          = AIF_part + Kep_Filter_Part;

end
