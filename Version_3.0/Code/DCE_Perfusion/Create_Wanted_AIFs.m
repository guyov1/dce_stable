function [ AIF_high_res, AIF_delayed_high_res, AIF, AIF_delayed] = Create_Wanted_AIFs( additonal_AIF_delay, time_vec_minutes_high_res,A1,sig1,T1,A2,sig2,T2,alpha,beta,s,tau, r_factor, Upsamp_factor )

AIF_high_res         = AIF_Parker(time_vec_minutes_high_res,A1,sig1,T1,A2,sig2,T2,alpha,beta,s,tau,0); %[mM]
AIF_delayed_high_res = AIF_Parker(time_vec_minutes_high_res,A1,sig1,T1,A2,sig2,T2,alpha,beta,s,tau, additonal_AIF_delay); %[mM]

% Downsample to wanted resolution
AIF                  = downsample(AIF_high_res, Upsamp_factor); %[mM]
AIF_delayed          = downsample(AIF_delayed_high_res, Upsamp_factor); %[mM]

% Normalize to 0->1 and Multiply by relaxivity factor to fit T1 values
AIF_high_res         = ReScale_AIF(r_factor, AIF_high_res);
AIF_delayed_high_res = ReScale_AIF(r_factor, AIF_delayed_high_res);
AIF                  = ReScale_AIF(r_factor, AIF);
AIF_delayed          = ReScale_AIF(r_factor, AIF_delayed);

end

