PKM3DFN=[WorkingP 'PKM3D.mat'];
load(PKM3DFN,'OutAIFParam','AMIdxs','ACXs','RMSs');
%%
TimeBetweenDCEVolsMin=TimeBetweenDCEVols/60;
InterpolationFactor=ceil(TimeBetweenDCEVols);
SampleTs=((1:nSVols)-1)*TimeBetweenDCEVolsMin;

HInterpolationFactor=ceil(InterpolationFactor*2);
Hdt=TimeBetweenDCEVolsMin/HInterpolationFactor;
HSampleTs=0:Hdt:SampleTs(end);

HConvIdxM=CreateConvIdxMFromSampleTs(numel(HSampleTs));
HTriB=HConvIdxM>0;
HConvIdxMTriB=HConvIdxM(HTriB);

nKeps=100;
Keps=gpowspace(0,15,nKeps,5)';
nKeps1=100;
Keps1I=floor(linspace(1,nKeps,nKeps1));

T1=0.17046;A1=0.809;sig1=0.0563;sig2=0.132;A2=0.330;T2=0.365;alpha=1.050;beta=0.1685;s=38.078;tau=0.483;
tauDelta=tau-T1;T2Delta=T2-T1;
% AIFParamsA=[T1 A1 sig1];
T1=1;

% AIF_Parker10t=@(x,t) AIF_Parker( t,A1,x(3)*sig1,x(1)*T1,A2*x(5),x(3)*sig2*x(6),x(1)*T1+T2Delta*x(4),alpha*x(7),beta*x(8),s*x(9),x(1)*T1+tauDelta*x(10) )*x(2)/1000;
% AIF_Parker8t=@(x,t) AIF_Parker( t,A1,x(3)*sig1,x(1)*T1,A2*x(5),x(3)*sig2*x(6),x(1)*T1+T2Delta*x(4),alpha*x(7),beta*x(8),s*1,x(1)*T1+tauDelta*0 )*x(2)/1000;
AIF_Parker8t=@(x,t) AIF_Parkerg2( t,A1,x(3)*sig1,x(1),A2*x(5),x(3)*sig2*x(6),x(1)+T2Delta*x(4),x(7),x(8))*x(2)/1000;

HAIF=AIF_Parker8t(OutAIFParam,HSampleTs);
tic;
HConvd=DCECostFuncgrT1ForConv(HAIF',Keps,HSampleTs,HConvIdxMTriB,HTriB);
tConv=toc;
disp(['Convolving the AIF took ' num2str(tConv)]);

ThreeSec=ceil(3/(Hdt*60));
TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;
nTDif=numel(TDif);

tic
SHConvd=zeros([nTDif nKeps numel(SampleTs)]);
HSHConvd=zeros([nTDif nKeps numel(HSampleTs)]);
SAIF=zeros([nTDif numel(SampleTs)]);
HSAIF=zeros([nTDif numel(HSampleTs)]);
for i=1:nTDif
    SHConvd(i,:,:)=interp1(HSampleTs,HConvd',SampleTs+TDif(i),[],'extrap')';
    SAIF(i,:)=interp1(HSampleTs,HAIF,SampleTs+TDif(i),[],'extrap');
    HSHConvd(i,:,:)=interp1(HSampleTs,HConvd',HSampleTs+TDif(i),[],'extrap')';
    HSAIF(i,:)=interp1(HSampleTs,HAIF,HSampleTs+TDif(i),[],'extrap');
end
tSampleConv=toc;
disp(['Sampling convolved CTCs took ' num2str(tSampleConv)]);
%%
MeanVol=loadniidata(MeanFN);

load([WorkingP 'CTCMsk.mat'],'CTC2DBigGood','MskCTCGood');
CTCFN=[WorkingP 'AfterCTC.mat'];
load(CTCFN,'DBrainMask');

MskCTCGood3D=DBrainMask;
MskCTCGood3D(MskCTCGood3D)=MskCTCGood;

Kep3D=MskCTCGood3D*0;
Kep3D(MskCTCGood3D)=Keps(AMIdxs(:,2));
BAT3D=MskCTCGood3D*0;
BAT3D(MskCTCGood3D)=-60*TDif(AMIdxs(:,1));
Vp3D=MskCTCGood3D*0;
Vp3D(MskCTCGood3D)=ACXs(1,:);
Ktrans3D=MskCTCGood3D*0;
Ktrans3D(MskCTCGood3D)=ACXs(2,:);
Ve3D=Ktrans3D./Kep3D;
Ve3D(~MskCTCGood3D)=0;
RMS3D=MskCTCGood3D*0;
RMS3D(MskCTCGood3D)=RMSs*1000;
% 
load(CTCFN,'CT1','B1');

Idx3DBig=NaN(size(Msk));
Idx3DBig(MskCTCGood3D)=1:sumn(MskCTCGood3D);

PrepareFN=[WorkingP 'AfterPrepare4D.mat'];
load(PrepareFN,'MaxRatioToBaseline');

MaxConcentration3D=NaN(size(Msk));
MaxConcentration3D(MskCTCGood3D)=max(CTC2DBigGood,[],2)*1000;

rRMS3D=RMS3D*NaN;
rRMS3D(MskCTCGood3D)=RMSs./max(CTC2DBigGood,[],2)';

Noise2D=EstimateNoise(CTC2DBigGood);
Noise3D=RMS3D*NaN;
Noise3D(MskCTCGood3D)=Noise2D;

RMStoNoise3D=RMS3D./(MaxConcentration3D.*max(Noise3D,1e-1));

Titles={'MeanVol' 'Kep','BAT','Vp','Ktrans','Ve' 'Msk' 'maxC' 'Noise' 'Enhancement%' 'RMS' 'rRMS3D' 'RMStoNoise' 'T1','B1'};
DCEResultsGUI({MeanVol Kep3D BAT3D Vp3D Ktrans3D Ve3D MskCTCGood3D MaxConcentration3D Noise3D (MaxRatioToBaseline-1)*100 RMS3D rRMS3D RMStoNoise3D CT1 B1},Titles,MskCTCGood3D,CTC2DBigGood,Idx3DBig,AMIdxs,ACXs,HSAIF,HSHConvd,Keps1I,SampleTs,HSampleTs,HAIF);

% Titles={'MeanVol' 'Kep','BAT','Vp','Ktrans','Ve' 'Msk'};
% DCEResultsGUI({MeanVol Kep3D BAT3D Vp3D Ktrans3D Ve3D Msk2
% },Titles,Msk2,CTC2D,Idx3D,AMIdxs,ACXs,HSAIF,HSHConvd,Keps1I,SampleTs,HSam
% pleTs,HAIF);