function Cost=AIFCostFuncM(WorkingP,AIFParams,DataToFit,SampleTs,HSampleTs,TDif,WPerCTC)

nTDif=numel(TDif);
AIF_Parker8t=@(x,t) AIF_Parkerg2( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8))*x(2);
HAIF=AIF_Parker8t(AIFParams,HSampleTs);
CHAIF=cumtrapz(HSampleTs,HAIF);
SAIF=zeros([nTDif numel(SampleTs)]);
CSAIF=zeros([nTDif numel(SampleTs)]);
for i=1:nTDif
    SAIF(i,:)=interp1(HSampleTs,HAIF,SampleTs+TDif(i),[],'extrap');
    CSAIF(i,:)=interp1(HSampleTs,CHAIF,SampleTs+TDif(i),[],'extrap');
end
[PKs Sims Cost] = FindPKBATgAIFMurase(WorkingP,DataToFit,SAIF,SampleTs,CSAIF,WPerCTC);