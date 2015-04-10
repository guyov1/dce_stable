
function to_fit = Tofts_model(time_vec_min, min_interval, AIF, Ktrans, Kep ,Vp)
   
    AIF_part        = Vp*AIF;
    
    if Sim_Struct.ignore_time_delta
        Kep_Filter_Part = filter(Ktrans*exp(-Kep*time_vec_min),1,AIF);
    else
        Kep_Filter_Part = filter(Ktrans*exp(-Kep*time_vec_min)*min_interval,1,AIF);
    end
    
    to_fit          = AIF_part + Kep_Filter_Part;

end
