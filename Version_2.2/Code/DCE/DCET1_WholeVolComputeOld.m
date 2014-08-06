[CVI BinCVI Bin2CVI]=ChooseRepVoxelsForAIFFind(Idx3D>0,CTC2D,BolusStart,5,5,2);

Options=struct('SubSecRes',{[]},'MaxTDif',{[]});

DataNoise=rmadCTC2D(CVI);
DataToFit=CTC2D(CVI,:);
PKMFN=[WorkingP 'PKM.mat'];

%%
load(PKMFN,'OutAIFParam');
%% Full 3D calculation
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
SAIF=zeros([nTDif numel(SampleTs)]);
for i=1:nTDif
    SHConvd(i,:,:)=interp1(HSampleTs,HConvd',SampleTs+TDif(i),[],'extrap')';
    SAIF(i,:)=interp1(HSampleTs,HAIF,SampleTs+TDif(i),[],'extrap');
end
tSampleConv=toc;
disp(['Sampling convolved CTCs took ' num2str(tSampleConv)]);
%%
clear R1F DCE2D2 adCTC2D B1 Baseline CT1 CosFA DCEM0

% load([WorkingP 'CTCBrainMsk.mat'],'CTC2DBrainMsk');
load([WorkingP 'CTCMsk.mat'],'CTC2DBigGood','MskCTCGood');
CTCFN=[WorkingP 'AfterCTC.mat'];
load(CTCFN,'DBrainMask');

MskCTCGood3D=DBrainMask;
MskCTCGood3D(MskCTCGood3D)=MskCTCGood;

N=size(CTC2DBigGood,1);
AMIdxs=zeros(N,2);
ACXs=zeros(2,N);
NAtATime=5000;
disp(['There are ' num2str(N) ' voxels to compute']);
before=now;
for i=1:NAtATime:N
    tic
    CurIs=i:min(N,i+NAtATime-1);
    [AMIdxs(CurIs,:) ACXs(:,CurIs) Sims] = FindKepBATgAIF(CTC2DBigGood(CurIs,:),SAIF,SHConvd,Keps1I);
    RMSs(CurIs)=  sqrt(mean((CTC2DBigGood(CurIs,:)-Sims).^2,2));
    t=toc;
    TimeFromStart=now-before;
    WillFinishAt=before+TimeFromStart*N/CurIs(end);
    disp(['Calculating ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
end
PKM3DFN=[WorkingP 'PKM3D.mat'];
save(PKM3DFN,'OutAIFParam','AMIdxs','ACXs','RMSs');
%% Output to Nii
MeanVol=loadniidata(MeanFN);

RMS3D=MskCTCGood3D*0;
RMS3D(MskCTCGood3D)=RMSs;

Kep3D=MskCTCGood3D*0;
Kep3D(MskCTCGood3D)=Keps(AMIdxs(:,2));
BAT3D=MskCTCGood3D*0;
BAT3D(MskCTCGood3D)=TDif(AMIdxs(:,1));
Vp3D=MskCTCGood3D*0;
Vp3D(MskCTCGood3D)=ACXs(1,:);
Ktrans3D=MskCTCGood3D*0;
Ktrans3D(MskCTCGood3D)=ACXs(2,:);
Ve3D=Ktrans3D./Kep3D;
Ve3D(~MskCTCGood3D)=0;

Raw2Nii(Kep3D,[WorkingP 'Kep.nii'],'float32', MeanFN);
Raw2Nii(BAT3D,[WorkingP 'BAT.nii'],'float32', MeanFN);
Raw2Nii(Vp3D,[WorkingP 'Vp.nii'],'float32', MeanFN);
Raw2Nii(Ktrans3D,[WorkingP 'Ktrans.nii'],'float32', MeanFN);
Raw2Nii(Ve3D,[WorkingP 'Ve.nii'],'float32', MeanFN);
Raw2Nii(RMS3D,[WorkingP 'RMS.nii'],'float32', MeanFN);