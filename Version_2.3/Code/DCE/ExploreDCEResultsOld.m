PKM3DFN=[WorkingP 'PKM3D' USStr '.mat'];
load(PKM3DFN,'OutAIFParam','PKs');
disp(PKM3DFN);
disp(OutAIFParam);
%%
TimeBetweenDCEVolsMin=TimeBetweenDCEVols/60;
InterpolationFactor=ceil(TimeBetweenDCEVols);
SampleTs=((1:nSVols)-1)*TimeBetweenDCEVolsMin;

HInterpolationFactor=ceil(InterpolationFactor*2);
Hdt=TimeBetweenDCEVolsMin/HInterpolationFactor;
HSampleTs=0:Hdt:SampleTs(end);

ThreeSec=ceil(3/(Hdt*60))*2;
TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;
nTDif=numel(TDif);

AIF_Parker8t=@(x,t) AIF_Parkerg2( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8))*x(2);

HAIF=AIF_Parker8t(OutAIFParam,HSampleTs);

%%
MeanVol=loadniidata(MeanFN);

load([WorkingP 'CTCMsk' USStr '.mat'],'CTC2DBigGood','MskCTCGood');
load(CTCFN,'DBrainMask');

MskCTCGood3D=DBrainMask;
MskCTCGood3D(MskCTCGood3D)=MskCTCGood;

Kep3D=MskCTCGood3D*0;
Kep3D(MskCTCGood3D)=PKs(:,2);
BAT3D=MskCTCGood3D*0;
BAT3D(MskCTCGood3D)=-60*TDif(PKs(:,1));
Vp3D=MskCTCGood3D*0;
Vp3D(MskCTCGood3D)=PKs(:,3);
Ktrans3D=MskCTCGood3D*0;
Ktrans3D(MskCTCGood3D)=PKs(:,4);
Ve3D=Ktrans3D./Kep3D;
Ve3D(~MskCTCGood3D)=0;
RMS3D=MskCTCGood3D*0;
RMS3D(MskCTCGood3D)=PKs(:,6)*1000;
F3D=MskCTCGood3D*0;
F3D(MskCTCGood3D)=max(0,PKs(:,7));
F03D=MskCTCGood3D*0;
F03D(MskCTCGood3D)=max(0,PKs(:,8));
RMSs=PKs(:,6)';
TP3D=MskCTCGood3D*0;
TP3D(MskCTCGood3D)=max(0,PKs(:,9));
TP3D(TP3D==2)=-0.5;
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
PKs(:,1)=-60*TDif(PKs(:,1));
Titles={'MeanVol' 'Kep','BAT','Vp','Ktrans','Ve' 'F' 'F0' 'TimeP' 'Msk' 'maxC' 'Noise' 'Enhancement%' 'RMS' 'rRMS3D' 'RMStoNoise' 'T1','B1'};
TTL=PKM3DFN;
DCEResultsGUI({MeanVol Kep3D BAT3D Vp3D Ktrans3D Ve3D -log10(F3D) -log10(F03D) TP3D MskCTCGood3D MaxConcentration3D Noise3D (MaxRatioToBaseline-1)*100 RMS3D rRMS3D RMStoNoise3D CT1 B1},Titles,MskCTCGood3D,CTC2DBigGood,Idx3DBig,1,PKs,1,1,1,SampleTs,HSampleTs,HAIF,TTL);
% DCEResultsGUI({MeanVol Kep3D BAT3D Vp3D Ktrans3D Ve3D F3D MskCTCGood3D MaxConcentration3D Noise3D (MaxRatioToBaseline-1)*100 RMS3D rRMS3D RMStoNoise3D CT1 B1},Titles,MskCTCGood3D,CTC2DBigGood,Idx3DBig,1,PKs,HSAIF,HSHConvd,Keps1I,SampleTs,HSampleTs,HAIF);

% Titles={'MeanVol' 'Kep','BAT','Vp','Ktrans','Ve' 'Msk'};
% DCEResultsGUI({MeanVol Kep3D BAT3D Vp3D Ktrans3D Ve3D Msk2
% },Titles,Msk2,CTC2D,Idx3D,AMIdxs,ACXs,HSAIF,HSHConvd,Keps1I,SampleTs,HSam
% pleTs,HAIF);