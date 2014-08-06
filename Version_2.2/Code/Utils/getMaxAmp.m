function MaxAmp=getMaxAmp(CTC2D)
Mx=max(CTC2D,[],2);
SMx=sort(Mx);
MaxAmp=SMx(numel(Mx)-10);