function Sims=SimulateCTC(SampleTs,HSampleTs,PKs,HAIF)
HConvIdxM=CreateConvIdxMFromSampleTs(numel(HSampleTs));
HTriB=HConvIdxM>0;
HConvIdxMTriB=HConvIdxM(HTriB);

HConvd2=DCECostFuncgrT1ForConv(HAIF',PKs(:,2),HSampleTs,HConvIdxMTriB,HTriB);

Hdt=HSampleTs(2)-HSampleTs(1);
dt=diff(SampleTs(1:2));
% ThreeSec=ceil(3/(Hdt*60));
% TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;
nToFit=size(PKs,1);
for i=1:nToFit
    SHConvd2(i,:)=interp1(HSampleTs,HConvd2(i,:)',HSampleTs-PKs(i,1),[],'extrap')';
    SHSAIF(i,:)=interp1(HSampleTs,HAIF,HSampleTs-PKs(i,1),[],'extrap');
end
%

for i=1:nToFit
%     Regressors=[SAIF(PKs(c,1),:); squeeze(SHConvd2(PKs(c,1),c,:))'];
    Regressors=[SHSAIF(i,:); squeeze(SHConvd2(i,:))];
        % BAT Kep Vp Ktrans Ve RMS fVal
    Sims(i,:)=((Regressors')*(PKs(i,[3 4])'));        
end
