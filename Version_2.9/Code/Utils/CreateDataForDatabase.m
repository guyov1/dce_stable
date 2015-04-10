StudyName='RoAs_20080122';
% StudyName='WhSa_20070813';
BasePath=['/media/OS/DATA/' StudyName '/'];
WorkingP=BasePath;

SDCEFA=25;
SDCETR=5.6;
SSNR=15;
SM0=10000;

aa=load([BasePath 'AfterCTC.mat']);
CTC2D=aa.CTC2D;
FinalT1=aa.FinalT1;
Msk=aa.Msk2;
b=load([BasePath 'PKM3D.mat']);
load([BasePath 'PKM' '.mat']);
B4=aa.B3;
B4(B4)=~aa.ImagB;
PKs=b.PKs(B4,:);
%%
USStr='';
CTCFN=[WorkingP 'AfterCTC' USStr '.mat'];
disp('Loading..');
load(CTCFN);
DCET1_PKf;
%%
Hdt=0.5/60;
TimeVec=SampleTs/60;
AIF_Parker9t=@(x,t) AIF_Parkerg3( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8),x(9),max(TimeVec))*x(2);

InspectedAIFParamsFN=[BasePath 'InspectedAIFParams.mat'];
if(exist(InspectedAIFParamsFN,'file'))
    Tmp=load(InspectedAIFParamsFN);
    OutAIFParam=Tmp.InspectedParams;
end

SampleTs=(aa.SampleTsNoBefore)/60;
HSampleTs=(0:Hdt:SampleTs(end));

HAIF=AIF_Parker9t(OutAIFParam,HSampleTs);

HConvIdxM=CreateConvIdxMFromSampleTs(numel(HSampleTs));
HTriB=HConvIdxM>0;
HConvIdxMTriB=HConvIdxM(HTriB);

% BATfinal VpFinal KtransFinal Kepfinal VeFinal
BATIdx=1;
VpIdx=2;
KtransIdx=3;
KepIdx=4;

ThreeSec = ceil(Options.MaxTDif_ForWholeVOI/(Hdt*60));
TDif     = -Hdt*ThreeSec:Hdt:Hdt*ThreeSec;
nTDif    = numel(TDif);
%%
% Idxs=1:10;
GoodIdxsB=isfinite(PKs(:,BATIdx)) & PKs(:,BATIdx)>1 & PKs(:,BATIdx)<nTDif;
GoodIdxsF=find(GoodIdxsB);

MskC=aa.Msk2;
MskC(MskC)=GoodIdxsB;
save([BasePath 'MskC.mat'],'MskC');

NPerIter=5000;
StartIdxs=1:NPerIter:numel(GoodIdxsF);

HHSampleTs=0:Hdt:floor(HSampleTs(end));
HHAIF=interp1(HSampleTs,HAIF',HHSampleTs);
HHConvIdxM=CreateConvIdxMFromSampleTs(numel(HHSampleTs));
HHTriB=HHConvIdxM>0;
HHConvIdxMTriB=HHConvIdxM(HHTriB);

Sims=zeros(numel(GoodIdxsF),numel(HHSampleTs));
%%
for ii=1:numel(StartIdxs)
    disp([ii numel(StartIdxs)]);
    disp(datestr(now));
    IdxsG=StartIdxs(ii):(min(StartIdxs(ii)+NPerIter-1,numel(GoodIdxsF)));
    Idxs=GoodIdxsF(IdxsG);
    CurKeps=PKs(Idxs,KepIdx);
    CurKeps(isnan(CurKeps))=0;
    
    HHConvd2=DCECostFuncgrT1ForConv(HHAIF',CurKeps,HHSampleTs,HHConvIdxMTriB,HHTriB);
    for i=1:size(HHConvd2,1)
        HConvd2(i,:)=interp1(HHSampleTs,HHConvd2(i,:),HHSampleTs,[],'extrap');
    end
    
    CurBATs=PKs(Idxs,BATIdx);
    CurBATs(isnan(CurBATs))=1;
    CurBATs=TDif(CurBATs);
    for i=1:numel(Idxs)
        SHConvd2(i,:)=interp1(HHSampleTs,HConvd2(i,:)',HHSampleTs+CurBATs(i),[],'extrap')';
        SHSAIF(i,:)=interp1(HSampleTs,HAIF,HHSampleTs+CurBATs(i),[],'extrap');
    end
    %
    
    CurKtranses=PKs(Idxs,KtransIdx);
    CurKtranses(isnan(CurKtranses))=0;
    for i=1:numel(Idxs)
        %     Regressors=[SAIF(PKs(c,1),:); squeeze(SHConvd2(PKs(c,1),c,:))'];
        Regressors=[SHSAIF(i,:); squeeze(SHConvd2(i,:))];
        % BATfinal VpFinal KtransFinal Kepfinal VeFinal
        Sims(IdxsG(i),:)=((Regressors')*([PKs(Idxs(i),[VpIdx]) CurKtranses(i)]'));
    end
end
disp('Finished simulation');
save([BasePath 'SimCTCClean.mat'],'Sims');
%% To Signal and add noise
GeneralDataFN=[WorkingP 'Params.mat'];
load(GeneralDataFN);
% Sims=load([BasePath 'SimCTCClean.mat'],'Sims');
% Sims=Sims.Sims;
% To T1
Sims=1./(Sims*AIFAmpCoeff+repmat(1./FinalT1(MskC),1,size(Sims,2)));
% To signal Out=SPGRfM(T1s,PDs,FAs,TR)
for i=1:size(Sims,2)
    Sims(:,i)=SPGRfM(Sims(:,i),SM0,SDCEFA,SDCETR);
end
% Add noise
Sims=Sims.*(1+randn(size(Sims))/SSNR);
save([BasePath 'SimSigNoised.mat'],'Sims');
disp('Finished noising');