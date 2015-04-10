% WorkingP='\\fmri-t9\users\Moran\DCE\DCE_Duration\001_HAVLIN_HAIM_MORDECHAY\Study20140831_120511_day224_T4\DCE_38min\HaMo_20140831\';
% AIFP='\\fmri-t9\users\Moran\DCE_Duration\01_HAVLIN_HAIM_MORDECHAY_3\Study20141019_100533\DCE\DCE6-3\HaMo_20141019\';
AIFP='\\fmri-t9\users\Moran\DCE_Duration\ParamsForParkerAIF\';
%%
ExportMatToCsv;
%% Load basic stuff
load('Export.mat');
ImportMatFromCsv;
handles=Export;

GeneralDataFN=[UntilDelimiter(Export.TTL,filesep,-1) filesep 'Params.mat'];
load(GeneralDataFN);

UseBAT=true;

Hdt=diff(handles.HSampleTs);
Hdt=Hdt(1);

handles.SampleTs=handles.SampleTs(~isnan(handles.SampleTs));
DTs=diff(handles.SampleTs);
nRegTimePoints=find(DTs>DTs(1)*2,1);
nExtraTPs=size(handles.SampleTs,2)-nRegTimePoints;
handles.HSampleTs=[0:Hdt:handles.SampleTs(nRegTimePoints) handles.SampleTs(nRegTimePoints+1:end)];
Idxs=Export.Idxs;
HConvIdxM=CreateConvIdxMFromSampleTs(numel(handles.HSampleTs));
HTriB=HConvIdxM>0;
HConvIdxMTriB=HConvIdxM(HTriB);

% BATfinal VpFinal KtransFinal Kepfinal VeFinal
BATIdx=1;
VpIdx=2;
KtransIdx=3;
KepIdx=4;
CurKeps=handles.CurPKs(:,KepIdx);
CurKeps(isnan(CurKeps))=0;

CurVps=handles.CurPKs(:,VpIdx)/AIFAmpCoeff;
CurVps(isnan(CurVps))=0;

nVoxels=numel(CurVps);

qq=load([AIFP 'InspectedAIFParams.mat']);
ww=load([AIFP 'InspectedAIFParamsTime.mat']);
handles.GoodTs=ww.InspectedParamsTimeSamples;
AIF_Parker9t=@(x,t) AIF_Parkerg3( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8),x(9),max(handles.GoodTs))*x(2);
AIFFunc=@(x,t) AIF_Parker9t(x,t);
handles.HAIF=AIFFunc(qq.InspectedParams,handles.HSampleTs);

A1=0.809;A2=0.330;T1=0.17046;T2=0.364;
sig1=0.055;sig2=0.134;alpha=1.064;beta=0.166;
s=37.772;tau=0.482;
% Before there were really slighltly different numbers, from other place?
% A~mM.min, T~min; sig~min, alpha~mM, beta,s~1/min, tau~min
% Time stamp (in minutes) for every temporal point
TimeBetweenDCEVols=diff(Export.SampleTs(1:2))*60;
nSVols=nRegTimePoints;
BolusStartSec=40;

TimeBetweenDCEVolsMin=TimeBetweenDCEVols/60;
SampleTs=((1:nSVols)-1)*TimeBetweenDCEVolsMin;
T1x=BolusStartSec/60;
% Population average c(t) according to Parker's article.
C=AIF_Parker(handles.HSampleTs,A1,sig1,T1x,A2,sig2,T2+T1x-T1,alpha,beta,s,tau+T1x-T1)/10;

% figure;plot(handles.HSampleTs,[C; handles.HAIF]);



HHSampleTs=0:Hdt:handles.HSampleTs(end);
HHAIF=interp1(handles.HSampleTs,handles.HAIF',HHSampleTs);
HHConvIdxM=CreateConvIdxMFromSampleTs(numel(HHSampleTs));
HHTriB=HHConvIdxM>0;
HHConvIdxMTriB=HHConvIdxM(HHTriB);

HHConvd2=DCECostFuncgrT1ForConv(HHAIF',CurKeps,HHSampleTs,HHConvIdxMTriB,HHTriB);
clear HConvd2
for i=1:size(HHConvd2,1)
    HConvd2(i,:)=interp1(HHSampleTs,HHConvd2(i,:),handles.HSampleTs,'linear','extrap');
end
% HConvd2=DCECostFuncgrT1ForConv(handles.HAIF',CurKeps,handles.HSampleTs,HConvIdxMTriB,HTriB);

dt=diff(handles.SampleTs(1:2));
% ThreeSec=ceil(3/(Hdt*60));
% TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;

if(~UseBAT)
    handles.nTDif=1;
    handles.TDif=0;
    CurBATs=CurBATs*0;
end

CurBATs=handles.CurPKs(:,BATIdx)/-60;
CurBATs(isnan(CurBATs))=1;
clear SHConvd2 SHSAIF
for i=1:nVoxels
    SHConvd2(i,:)=interp1(handles.HSampleTs,HConvd2(i,:)',handles.HSampleTs+CurBATs(i),'linear','extrap')';
    SHSAIF(i,:)=interp1(handles.HSampleTs,handles.HAIF,handles.HSampleTs+CurBATs(i),'linear','extrap');
end
%
CurKtranses=handles.CurPKs(:,KtransIdx)/AIFAmpCoeff;
CurKtranses(isnan(CurKtranses))=0;
clear Sims
for i=1:nVoxels
    Regressors=[SHSAIF(i,:); squeeze(SHConvd2(i,:))];
    Sims(i,:)=((Regressors')*([handles.CurPKs(i,[VpIdx]) CurKtranses(i)]'));        
end

%
SAIF=zeros([handles.nTDif numel(handles.SampleTs)]);
CSAIF=zeros([handles.nTDif numel(handles.SampleTs)]);
CHAIF = cumtrapz(handles.HSampleTs,handles.HAIF);
for i=1:handles.nTDif
    SAIF(i,:)=interp1(handles.HSampleTs,handles.HAIF,handles.SampleTs+handles.TDif(i),'linear','extrap');
    CSAIF(i,:)=interp1(handles.HSampleTs,CHAIF,handles.SampleTs+handles.TDif(i),'linear','extrap');
end
% Simulate
clear SSims
for i=1:size(Sims,1)
    SSims(i,:)=interp1(handles.HSampleTs,Sims(i,:),handles.SampleTs,'linear','extrap');
end
SSimsA=SSims;
%% Extract for simulated
SSimsB=SSims.*(1+randn(size(SSims))*0.2); %insert noise (0.2=20%)
LastTimePoints=handles.SampleTs(end-nExtraTPs:end);
ShortName=handles.TTL(end-22:end-10);
clear PKs_check
%    1        2         3        4          5     6        7    8   9    10
% BATfinal VpFinal KtransFinal Kepfinal VeFinal RSSFinal RSS0 RSS1 RSS2 RSS3
for k=1:nExtraTPs+1
    PKs_check(k,:,:) = FindPKBATgAIFMuraseF4Models_TProb(SSimsB(:,1:end-nExtraTPs-1+k),SAIF(:,1:end-nExtraTPs-1+k),handles.SampleTs(:,1:end-nExtraTPs-1+k),CSAIF(:,1:end-nExtraTPs-1+k));
end

clear SVps SKtranss SKeps
for k=1:nExtraTPs+1
    SVps(:,k)=PKs_check(k,:,VpIdx);
    SKtranss(:,k)=PKs_check(k,:,KtransIdx);
    %SKeps(:,k)=PKs_check(k,:,KepIdx); Run with Model Selection
    SKeps(:,k)=PKs_check(k,:,22); %Run without Model Selection
end
%%
figure;
subplot(1,3,1);
plot(LastTimePoints,SVps','*');hold on;plot(LastTimePoints,repmat(CurVps,[1 nExtraTPs+1])','-')
title('Real values-r, Estimated values-b')
xlabel('minutes');
ylabel('v_p')
subplot(1,3,2);
plot(LastTimePoints,SKtranss','*');hold on;plot(LastTimePoints,repmat(CurKtranses,[1 nExtraTPs+1])','-')
title(ShortName)
xlabel('minutes');
ylabel('K^{trans}')
subplot(1,3,3);
plot(LastTimePoints,SKeps','*');hold on;plot(LastTimePoints,repmat(CurKeps,[1 nExtraTPs+1])','-')
xlabel('minutes');
ylabel('k_{ep}')
%%
Stuff.HSampleTs=handles.HSampleTs;
Stuff.HAIF=handles.HAIF;
Stuff.HHAIF=HHAIF;
Stuff.HHSampleTs=HHSampleTs;
Stuff.HHConvIdxMTriB=HHConvIdxMTriB;
Stuff.HHTriB=HHTriB;
%%
clear ASims
for CurVox=1:size(PKs_check,2)
    ASims(:,:,CurVox)=PK2CTC(squeeze(PKs_check(:,CurVox,1:4)),Stuff);
end
%%
figure;
WhichToShow=1; %:size(PKs_check,2);
%plot(Stuff.HSampleTs,Stuff.HAIF/200,'k');hold on;
plot(Stuff.HSampleTs,Stuff.HAIF/10,'k');hold on;
for CurVox=WhichToShow
    plot(handles.HSampleTs,squeeze(ASims(:,:,CurVox))');hold on;
    plot(Export.SampleTs(1:size(Export.CurCTCs,2)),Export.CurCTCs(CurVox,:),'ko');
end
%% Export to csv
CurVs=[CurVps CurKtranses CurKeps];
EstM=[LastTimePoints;SVps;SKtranss;SKeps];
clear OutM;
OutM(1,1)=nVoxels;
OutM(2:(nVoxels+1),1:3)=CurVs;
OutM((nVoxels+2):(nVoxels*4+2),1:numel(LastTimePoints))=EstM;
csvwrite(['SimLong_' ShortName '.csv'],OutM);