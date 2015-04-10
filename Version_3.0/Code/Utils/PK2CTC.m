function Sims=PK2CTC(CurPKs,Stuff)


BATIdx=1;
VpIdx=2;
KtransIdx=3;
KepIdx=4;

CurKeps=CurPKs(:,KepIdx);
CurKeps(isnan(CurKeps))=0;

HHConvd2=DCECostFuncgrT1ForConv(Stuff.HHAIF',CurKeps,Stuff.HHSampleTs,Stuff.HHConvIdxMTriB,Stuff.HHTriB);
clear HConvd2
for i=1:size(HHConvd2,1)
    HConvd2(i,:)=interp1(Stuff.HHSampleTs,HHConvd2(i,:),Stuff.HSampleTs,'linear','extrap');
end

CurBATs=CurPKs(:,BATIdx)/-60;
CurBATs(isnan(CurBATs))=1;
clear SHConvd2 SHSAIF
for i=1:size(CurPKs,1)
    SHConvd2(i,:)=interp1(Stuff.HSampleTs,HConvd2(i,:)',Stuff.HSampleTs+CurBATs(i),'linear','extrap')';
    SHSAIF(i,:)=interp1(Stuff.HSampleTs,Stuff.HAIF,Stuff.HSampleTs+CurBATs(i),'linear','extrap');
end
%
CurKtranses=CurPKs(:,KtransIdx);
CurKtranses(isnan(CurKtranses))=0;
clear Sims
for i=1:size(CurPKs,1)
    Regressors=[SHSAIF(i,:); squeeze(SHConvd2(i,:))];
    Sims(i,:)=((Regressors')*([CurPKs(i,[VpIdx]) CurKtranses(i)]'));        
end