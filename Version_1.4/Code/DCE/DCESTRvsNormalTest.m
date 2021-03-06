WorkingP='/data/home/gilad/Desktop/DCE/gilad/Database/DCEOut/KaRo_20070911/';

a=load([WorkingP 'AfterCTC.mat']);

nVols=size(a.CTC2D,2);
TimeBetweenDCEVols=a.TimeBetweenDCEVols;
STimeBetweenDCEVols=a.TimeBetweenDCEVols;
STimeBetweenDCEVolsMin=a.TimeBetweenDCEVols/60;
SMin2Timestep=1/STimeBetweenDCEVolsMin;
nSVols=nVols;
SampleTs=((1:nSVols)-1)*STimeBetweenDCEVolsMin;

ManualArt3D=loadniidata([WorkingP 'manualArt.nii']);
F1=find(ManualArt3D);
F2=find(a.Msk2);
CVI=find(ismember(F2,F1));
ArtCTCs=a.CTC2D(CVI,:);
MaxAmp=max(ArtCTCs(:));
DataNoise=EstimateNoise(ArtCTCs);
[Tmp, BolusStart]=max(ArtCTCs(1,:));
%%
Sdt=0.01;
SEndTime=10; % minutes
Sts=-SEndTime:Sdt:SEndTime; %(-nSVols*2:dVols:nSVols*2)./Min2Timestep;
options = struct('GradObj','off','Display','off');
InterestingTimePoints=Sts>-1;

Sts=SampleTs;

% A~mM.min, T~min; sig~min, alpha~mM, beta,s~1/min, tau~min
T1=0.17046;A1=0.809;sig1=0.0563;sig2=0.132;A2=0.330;T2=0.365;alpha=1.050;beta=0.1685;s=38.078;tau=0.483;
tauDelta=tau-T1;T2Delta=T2-T1;
% AIFParamsA=[T1 A1 sig1];
T1=1;
% AIF_Parker3=@(x) AIF_Parker(Sts,A1,x(3)*sig1,x(1)*T1,A2,x(3)*sig2,x(1)*T1+T2Delta,alpha,beta,s,x(1)*T1+tauDelta)*x(2)/1000;
% AIF_Parker6=@(x) AIF_Parker(Sts,A1,x(3)*sig1,x(1)*T1,A2*x(5),x(3)*sig2*x(6),x(1)*T1+T2Delta*x(4),alpha,beta,s,x(1)*T1+tauDelta)*x(2)/1000;
% AIF_Parker10=@(x) AIF_Parker( Sts,A1,x(3)*sig1,x(1)*T1,A2*x(5),x(3)*sig2*x(6),x(1)*T1+T2Delta*x(4),alpha*x(7),beta*x(8),s*x(9),x(1)*T1+tauDelta*x(10) )*x(2)/1000;
% OrigParkerAIF=AIF_Parker3(ones(1,3));

ParamAIFCoeff=[0.8,25,3;... % First bolus time
    0.25,4,1;... % First bolus magnitude
    0.25,10,2;... % First bolus width
    0.5,10,1;... % 2nd bolus time
    0.25,4,1;... % 2nd bolus magnitude
    0.25,4,1;... % 2nd bolus width
    0.5,2,1;... % alpha
    0.5,2,1;... % beta
    0.5,2,1;... % s
    0.5,5,1]; % tauDelta

ParamAIFCoeff(1,3)=BolusStart*TimeBetweenDCEVols/60;
ParamAIFCoeff(3,:)=[1 5 3];
ParamAIFCoeff(1,1)=0.1;
ParamAIFCoeff(1,2)=SampleTs(end)/2;
ParamAIFCoeff(5,1)=0;
ParamAIFCoeff(7,:)=[0.2 2 0.4];
ParamAIFCoeff(8,:)=[0 0.3 beta];
% ParamAIFCoeff(7,2)=20;
% ParamAIFCoeff(8,2)=20;
ParamAIFCoeff(10,1)=0;
ParamAIFCoeff(10,3)=0.5;
ParamAIFCoeff(2,2)=10;
ParamAIFCoeff(6,:)=[1 5 3];

nAIFParams=size(ParamAIFCoeff,1);

Defaults.SubSampling=1;
Defaults.nVolsToRemoveFromEnd=0;
Defaults.SubSecondResolution=2;
Defaults.MinFirstBolusStd=2;
Defaults.EM_Num_Of_Iterations=5;
Defaults.FMS_TolFun=1e-11;
Defaults.FMS_MaxFunEvals=10000;
Defaults.FMS_MaxIter=10000;
Defaults.MaxTDif_ForAIFSearch=3;
Defaults.MaxTDif_ForWholeVOI=6;
Defaults.Rep_MaxAroundBolus=10;
Defaults.Rep_RatioToEnd=10;
Defaults.Rep_nPerSet=1;
Options=Defaults;
%%
USStr='';
AIFFinderFN=[WorkingP 'STRvsNoralAIFFindData' USStr '.mat'];
[PKOut OutAIFParam]=AIFTryf2(ArtCTCs,MaxAmp,DataNoise,ParamAIFCoeff,nSVols,TimeBetweenDCEVols,Defaults,false,AIFFinderFN);
%%
TimeBetweenDCEVolsMin=TimeBetweenDCEVols/60;
InterpolationFactor=ceil(TimeBetweenDCEVols);
SampleTs=((1:nSVols)-1)*TimeBetweenDCEVolsMin;

HInterpolationFactor=ceil(InterpolationFactor*Options.SubSecondResolution);
% HInterpolationFactor=ceil(InterpolationFactor*Options.SubSecRes);
Hdt=TimeBetweenDCEVolsMin/HInterpolationFactor;
HSampleTs=0:Hdt:SampleTs(end);

ThreeSec=ceil(Options.MaxTDif_ForWholeVOI/(Hdt*60));
TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;
nTDif=numel(TDif);

AIF_Parker8t=@(x,t) AIF_Parkerg2( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8))*x(2);

HAIF=AIF_Parker8t(OutAIFParam,HSampleTs);

CHAIF=cumtrapz(HSampleTs,HAIF);
    
SAIF=zeros([nTDif numel(SampleTs)]);
CSAIF=zeros([nTDif numel(SampleTs)]);
for i=1:nTDif
    SAIF(i,:)=interp1(HSampleTs,HAIF,SampleTs+TDif(i),[],'extrap');
    CSAIF(i,:)=interp1(HSampleTs,CHAIF,SampleTs+TDif(i),[],'extrap');
end
%%
figure;plot(SampleTs,NormalizeByRows(ArtCTCs)*MaxAmp,'*-');hold on;plot(HSampleTs,HAIF,'m-','LineWidth',2);
%% Create Sims
nToFit=size(ArtCTCs,1);
PKs=PKOut(:,[2 1 4 3]);

WorkP2='/data/home/gilad/Desktop/DCE/gilad/Database/DCEOut/KaRo_20071111/';
PKM3DFN=[WorkP2 'PKM3D' USStr '.mat'];
d=load(PKM3DFN,'PKs');
%%
SPKs=max(d.PKs.*(randn(size(d.PKs))*0.1+1),0);

MinVp=0;
MaxVp=1;
MinVe=0.05;
MaxVe=10; % Account for AIF maxAmp too big -> ktrans small
MaxKtrans=2.2;
MaxKep=MaxKtrans/MinVe;
% BAT Kep Vp Ktrans Ve RMS fVal f0Val
SPKs(:,1)=TDif(d.PKs(:,1))'+rand(size(SPKs,1),1)/60-(1/120);
SPKs(:,1)=(rand(size(SPKs,1),1)*6-3)/60;
SPKs(:,2)=min(MaxKep,SPKs(:,2));
SPKs(:,3)=min(MaxVp,SPKs(:,3));
SPKs(:,4)=min(MaxKtrans,SPKs(:,4));
%%
CVI=getKrandomSamples(size(SPKs,1),1000);
PKs=SPKs(CVI,:);
HSims=SimulateCTC(SampleTs,HSampleTs,PKs,HAIF);
% HSims=HSims.*(randn(size(HSims))*0.05+1);
HSims=HSims+randn(size(HSims))*0.00007;
Sims=interp1(HSampleTs,HSims',SampleTs,[],'extrap')';
figure;plot(HSampleTs,HSims'); hold on; plot(SampleTs,Sims,'.-');% plot(SampleTs,ArtCTCs,'*');
%% Using STR
Tmpstr= FindPKBATgAIFMuraseF4Models_TProb(Sims,SAIF,SampleTs,CSAIF);
%%
% STR_PKs=Tmpstr(:,[19 22 20 21]);
STR_PKs=Tmpstr(:,[1 4 2 3]);
STR_PKs(isfinite(STR_PKs(:,1)),1)=TDif(STR_PKs(isfinite(STR_PKs(:,1)),1));
figure(99);clf;
for i=1:4
    subplot(2,2,i);
    plot(PKs(:,i),STR_PKs(:,i)./PKs(:,i),'*');
    BB=STR_PKs(:,i);
%     BB(BB<0.01)=Inf;
    BB(PKs(:,i)<=0.1)=Inf;
    
%     title(gCost(PKs(isfinite(BB),i),BB(isfinite(BB)),'Corr'));
    title(gCost(PKs(isfinite(BB),i),BB(isfinite(BB)),{'RMS','Relative'}));
end
%% Using normal
figure;clf;hold on;
plot(HSampleTs,HAIF,'k','LineWidth',7);
for k=1:(HInterpolationFactor-1)
    plot(SampleTs+TimeToStart(k),SmpAIF(k,:),'Color',ClrM(k,:));
end
%%
for k=1:(HInterpolationFactor-1)
    CurSims=interp1(HSampleTs,HSims',SampleTs+TimeToStart(k)*1,[],'extrap')';

    SmpAIF(k,:)=interp1(HSampleTs,HAIF,SampleTs+TimeToStart(k),[],'extrap');
    % figure;plot(HSampleTs*60,HAIF,'b.');hold on;plot((SampleTs+TimeToStart(k))*60,SmpAIF(k,:),'g-*');
    SCSAIF=cumtrapz(SampleTs,SmpAIF(k,:));
    Tmp = FindPKBATgAIFMuraseF4Models_TProb(CurSims,SmpAIF(k,:),SampleTs,SCSAIF);
    AllTmps(k,:,:)=Tmp;
end
%%
clear SmpPKs
TimeToStart=linspace(0,dt,HInterpolationFactor);
ClrM=jet(HInterpolationFactor);
figure(199);clf
for k=1:(HInterpolationFactor-1)
    Tmp=squeeze(AllTmps(k,:,:));
%     SmpPKs(k,:,:)=Tmp(:,[19 22 20 21]);
    SmpPKs(k,:,:)=Tmp(:,[1 4 2 3]);
    %

    for i=2:4
        subplot(2,2,i);hold on;
        plot(PKs(:,i),squeeze(SmpPKs(k,:,i))'./PKs(:,i),'*','Color',ClrM(k,:));
        BB=squeeze(SmpPKs(k,:,i));
        BB(PKs(:,i)<=0.1)=Inf;
%         SmpCorrs(k,i)=gCost(PKs(isfinite(BB),i),squeeze(SmpPKs(k,isfinite(BB),i))','Corr');
        SmpCorrs(k,i)=gCost(PKs(isfinite(BB),i),squeeze(SmpPKs(k,isfinite(BB),i))',{'RMS','Relative'});
    end
end
subplot(2,2,1);
plot(SmpCorrs);