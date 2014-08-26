function rmadCTC2D=EstimateNoise(CTC2D)
adCTC2D=abs(diff(CTC2D,[],2));
madCTC2D=median(adCTC2D,2);
rmadCTC2D=madCTC2D./abs(max(CTC2D,[],2));
