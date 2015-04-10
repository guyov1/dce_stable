function [PKOut OutAIFParam]=AIFTryf(WorkingP,DataToFit,BinCVI,DataNoise,ParamAIFCoeff,nSVols,TimeBetweenDCEVols,Options,ShowFig,AIFFinderFN)
TimeBetweenDCEVolsMin=TimeBetweenDCEVols/60;
InterpolationFactor=ceil(TimeBetweenDCEVols);
SampleTs=((1:nSVols)-1)*TimeBetweenDCEVolsMin;

Defaults=struct('nIterations',10,'MaxTDif',3,'SubSecRes',2);
options = optimset('TolFun',1e-11,'MaxFunEvals',10000,'MaxIter',10000);

Options=ExtendStruct(Options,Defaults);
%%

nToFit=size(DataToFit,1);

T1=0.17046;A1=0.809;sig1=0.0563;sig2=0.132;A2=0.330;T2=0.365;alpha=1.050;beta=0.1685;s=38.078;tau=0.483;
tauDelta=tau-T1;T2Delta=T2-T1;
% AIFParamsA=[T1 A1 sig1];
T1=1;

% AIF_Parker10t=@(x,t) AIF_Parker( t,A1,x(3)*sig1,x(1)*T1,A2*x(5),x(3)*sig2*x(6),x(1)*T1+T2Delta*x(4),alpha*x(7),beta*x(8),s*x(9),x(1)*T1+tauDelta*x(10) )*x(2)/1000;
% AIF_Parker8t=@(x,t) AIF_Parker(
% t,A1,x(3)*sig1,x(1)*T1,A2*x(5),x(3)*sig2*x(6),x(1)*T1+T2Delta*x(4),alpha*x(7),beta*x(8),s*1,x(1)*T1+tauDelta*0 )*x(2)/1000;
AIF_Parker8t=@(x,t) AIF_Parkerg2( t,A1,x(3)*sig1,x(1),A2*x(5),x(3)*sig2*x(6),x(1)+T2Delta*x(4),x(7),x(8))*x(2)/1000;
HInterpolationFactor=ceil(InterpolationFactor*Options.SubSecRes);
Hdt=TimeBetweenDCEVolsMin/HInterpolationFactor;
HSampleTs=0:Hdt:SampleTs(end);
nKeps=100;
Keps=gpowspace(0,15,nKeps,5)';
nKeps1=100;
Keps1I=floor(linspace(1,nKeps,nKeps1));
WReg=10000;
WRegBolus=10000;
WSmooth1=100000;
WSmooth1Bolus=1000;
WSmooth2=100000;
WSmooth2Bolus=100;
% TDif=(-0*TimeBetweenDCEVolsMin):Hdt:(2*TimeBetweenDCEVolsMin);
ThreeSec=ceil(Options.MaxTDif/(Hdt*60));
TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;
nTDif=numel(TDif);
ZeroTIdx=find(TDif==0);
OldParams=ParamAIFCoeff(1:8,3);
% 
% ParamAIFCoeff=[0.8,25,3;... % First bolus time
%     0.25,4,1;... % First bolus magnitude
%     0.25,10,2;... % First bolus width
%     0.5,10,1;... % 2nd bolus time
%     0.25,4,1;... % 2nd bolus magnitude
%     0.25,4,1;... % 2nd bolus width
%     0.5,2,1;... % alpha
%     0.5,2,1;... % beta
%     0.5,2,1;... % s
%     0.5,5,1]; % tauDelta

ParamAIFCoeff(3,1)=1;
ParamAIFCoeff(2,[1 2])=[0.01 1000];
%% Start with AIF from Artery
VesselsIdxs=find(BinCVI==max(BinCVI));
MaxAmp=max(max(DataToFit(VesselsIdxs,:)));
NormalizedVessels=NormalizeByRows(DataToFit(VesselsIdxs,:))*MaxAmp;
MeanVessel=mean(NormalizedVessels,1);
AIF_Parker8tx=@(x,t) NormalizeByRows(AIF_Parker8t(x,t)).*MaxAmp;
OldParams = lsqcurvefit(AIF_Parker8tx,OldParams(1:8)',SampleTs,MeanVessel,ParamAIFCoeff(1:8,1), ParamAIFCoeff(1:8,2),options)';
HAIF=AIF_Parker8t(OldParams,HSampleTs);
OldParams(2)=OldParams(2)*MaxAmp/max(HAIF);
HAIF=AIF_Parker8t(OldParams,HSampleTs);
% figure(99);plot(SampleTs,NormalizedVessels,'k-*');hold on;plot(HSampleTs,HAIF,'b','LineWidth',3);
%%
tic
% HBolusStart=(BolusStart-1)*HInterpolationFactor+1;
HBolusStart=ceil(OldParams(1)/Hdt);
nTimePoints=numel(SampleTs);
nHTimePoints=numel(HSampleTs);
ConvIdxM=CreateConvIdxMFromSampleTs(numel(SampleTs));
TriB=ConvIdxM>0;
% ConvIdxMTriB=ConvIdxM(TriB);
HConvIdxM=CreateConvIdxMFromSampleTs(numel(HSampleTs));
HTriB=HConvIdxM>0;
HConvIdxMTriB=HConvIdxM(HTriB);

WPerCTC=1./(max(DataToFit,[],2).*DataNoise);
HAIF=AIF_Parker8t(OldParams,HSampleTs);
tBase=toc;
disp(['Preparation took ' num2str(tBase)]);
%%
IterCont=true;
IterCount=0;
clear AllAIF2s AllAIF2cs AllAIFPs AllRMSs ARMS
AllAIFPs(1,:)=OldParams;
while(IterCont)
    IterCount=IterCount+1;
    IterCont=IterCount<Options.nIterations;
    disp(['Iteration #' num2str(IterCount)]);
    %
    tic;
    HConvd=DCECostFuncgrT1ForConv(HAIF',Keps,HSampleTs,HConvIdxMTriB,HTriB);
    tConv=toc;
    disp(['Convolving the AIF took ' num2str(tConv)]);
%     if(ShowFig)
%         figure;plot(HSampleTs,NormalizeByRows(HConvd(:,:))')
%     end
    %
    tic
    SHConvd=zeros([nTDif nKeps numel(SampleTs)]);
    HSHConvd=zeros([nTDif nKeps numel(HSampleTs)]);
    SAIF=zeros([nTDif numel(SampleTs)]);
    HSAIF=zeros([nTDif numel(HSampleTs)]);
    CHAIF=cumsum(HAIF);
    CSAIF=zeros([nTDif numel(SampleTs)]);
    for i=1:nTDif
        SHConvd(i,:,:)=interp1(HSampleTs,HConvd',SampleTs+TDif(i),[],'extrap')';
        SAIF(i,:)=interp1(HSampleTs,HAIF,SampleTs+TDif(i),[],'extrap');
        CSAIF(i,:)=interp1(HSampleTs,CHAIF,SampleTs+TDif(i),[],'extrap');
        HSHConvd(i,:,:)=interp1(HSampleTs,HConvd',HSampleTs+TDif(i),[],'extrap')';
        HSAIF(i,:)=interp1(HSampleTs,HAIF,HSampleTs+TDif(i),[],'extrap');
    end
    tSampleConv=toc;
    disp(['Sampling convolved CTCs took ' num2str(tSampleConv)]);
%     if(ShowFig)
%         figure;plot(SampleTs,NormalizeByRows(squeeze(SHHConvd(:,50,:)))')
%     end
    % 
    tic
%     [MIdxs CXs Sims] = FindKepBATgAIF(DataToFit,SAIF,SHConvd,Keps1I);
    CumSumAIFdTSamp=CSAIF*TimeBetweenDCEVolsMin;
    [PKs Sims] = FindPKBATgAIFMurase(WorkingP,DataToFit,SAIF,SampleTs,CumSumAIFdTSamp);
    
    tFindPK=toc;
    disp(['Finding Kep, BAT took ' num2str(tFindPK)]);
    ARMS(IterCount)=sqrt(mean(((Sims-DataToFit).^2),2))'*WPerCTC;
    % Deconv
    tic
    RHS=[];
    LHS=[];

    [U, QQQ, IB]=unique(MIdxs(:,2));
    nUKep=numel(U);
    EKepM=zeros(nHTimePoints);
    EKepMs=zeros([nUKep nTimePoints nHTimePoints]);
    for i=1:nUKep
        EKepV=exp(-Keps(Keps1I(U(i)))*HSampleTs)*Hdt;
        EKepM(HTriB)=EKepV(HConvIdxMTriB);
        EKepMs(i,:,:)=EKepM(1:HInterpolationFactor:end,:); %
    end
    TShift=MIdxs(:,1)-ZeroTIdx;
    LittleI=eye(nHTimePoints-min(0,min(TShift)));
    LittleIBase=LittleI(1:HInterpolationFactor:end,:);
    LHS=[];
    for c=1:nToFit
        EKepM=squeeze(EKepMs(IB(c),:,:));
        LittleI=LittleIBase(1:nTimePoints,1:nHTimePoints);
        if(TShift(c)>=0)
            EKepM(:,(TShift(c)+1):end)=EKepM(:,1:end-TShift(c));
            EKepM(:,1:TShift(c))=0;
            LittleI(:,(TShift(c)+1):end)=LittleI(:,1:end-TShift(c));
            LittleI(:,1:TShift(c))=0;
        else
            EKepM(:,1:end+TShift(c))=EKepM(:,(1-TShift(c)):end);
            EKepM(:,end+TShift(c)+1:end)=0;
            LittleI(:,1:end+TShift(c))=LittleI(:,(1-TShift(c)):end);
            LittleI(:,end+TShift(c)+1:end)=0;
        end

        CurLHS=EKepM*CXs(2,c)+LittleI*CXs(1,c);
        GoodRows=sum(LittleI,2)>0;
        CurLHS=CurLHS*WPerCTC(c);

        LHS=[LHS; CurLHS(GoodRows,:)];
        CurRHS=DataToFit(c,:)'*WPerCTC(c);

        RHS=[RHS; CurRHS(GoodRows,:)];
    end
    % Regularization
    AroundBolus=(1:nHTimePoints);
    AroundBolus=AroundBolus<HBolusStart+HInterpolationFactor*1 & AroundBolus>HBolusStart-HInterpolationFactor*1;
    LHSReg=diag(AroundBolus*WRegBolus+(~AroundBolus)*WReg);
    RHSReg=HAIF'*WReg;
    LHSSmooth1=diag(AroundBolus*WSmooth1Bolus+(~AroundBolus)*WSmooth1);
    LHSSmooth1=LHSSmooth1(1:end-1,:)-LHSSmooth1(2:end,:);
    RHSSmooth1=LHSSmooth1*HAIF';%  zeros(nHTimePoints-1,1);
    LHSSmooth2=diag(AroundBolus*WSmooth2Bolus+(~AroundBolus)*WSmooth2);
    LHSSmooth2=LHSSmooth2(1:end-2,:)-2*LHSSmooth2(2:end-1,:)+LHSSmooth2(3:end,:);
    RHSSmooth2=LHSSmooth2*HAIF';
    LHS=[LHS; LHSReg; LHSSmooth1; LHSSmooth2];
    RHS=[RHS; RHSReg; RHSSmooth1; RHSSmooth2];
    %
    AIF2=LHS\RHS;
%     AIF2=lsqnonneg(LHS,RHS,AIF2);
%     AIF2=blocknnls(LHS,RHS);
    RMS1=gCost(LHS*AIF2,RHS,'RMS');
    tLS=toc;
    disp(['Finding AIF by LS took ' num2str(tLS) ' to RMS ' num2str(RMS1)]);
    %
    tic
%     AIF_Func1=@(x) AIF_Parker10t(x,HSampleTs);
    AIF_Func1=@(x) AIF_Parker8t(x,HSampleTs);
%     NewPc=FindAIFCoeffs(AIF2', AIF_Func1, OldParams, ParamAIFCoeff(:,1), ParamAIFCoeff(:,2))';
%     NewPc =
%     lsqcurvefit(AIF_Parker10t,OldParams,HSampleTs,AIF2',ParamAIFCoeff(:,1), ParamAIFCoeff(:,2),options)';
    NewPc = lsqcurvefit(AIF_Parker8t,OldParams(1:8),HSampleTs,AIF2',ParamAIFCoeff(1:8,1), ParamAIFCoeff(1:8,2),options)';
%     x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options)
    NewPc(2)=OldParams(2);
    NewPc(1)=OldParams(1);
    AIF2c=AIF_Func1(NewPc)';
    RMS2=gCost(LHS*AIF2c,RHS,'RMS');
    tToFuncParams=toc;
    disp(['Converting to functional form took ' num2str(tToFuncParams) ' to RMS ' num2str(RMS2)]);
    if(ShowFig)
        figure;plot(HSampleTs,[HAIF' AIF2 AIF2c])
    end
    OldParams=NewPc;
    AllAIF2s(IterCount,:)=AIF2;
    AllAIF2cs(IterCount,:)=AIF2c;
    AllAIFPs(IterCount+1,:)=NewPc;
    AllRMSs(IterCount,:)=[RMS1 RMS2];
%     HAIF=AIF_Parker10t(OldParams,HSampleTs);
    HAIF=AIF_Parker8t(OldParams,HSampleTs);
    HBolusStart=ceil(OldParams(1)/Hdt);
%     HAIF=AIF2;
end
disp('Finished iterations');
if(nargin>8)
    save(AIFFinderFN,'AllAIF2s','AllAIF2cs','AllRMSs','ARMS','AllAIFPs');
end
%%
[Tmp, MinARMSI]=min(ARMS)
% MinARMSI=1;
%%
if(ShowFig)
    figure(7234723);clf;
    plot(AllAIF2s');

    figure(7234724);clf;
    plot(AllAIF2cs');hold on;
    plot(AllAIF2cs(MinARMSI-1,:),'LineWidth',2);

    figure(9282);clf;
    plot([AllRMSs ARMS'/30]);
end
%%
OutAIFParam=AllAIFPs(MinARMSI,:);
HAIF=AIF_Parker8t(AllAIFPs(MinARMSI,:),HSampleTs);
tic;
HConvd=DCECostFuncgrT1ForConv(HAIF',Keps,HSampleTs,HConvIdxMTriB,HTriB);
tConv=toc;
disp(['Convolving the AIF took ' num2str(tConv)]);
if(ShowFig)
    figure;plot(HSampleTs,NormalizeByRows(HConvd(:,:))')
end
%
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
if(ShowFig)
    figure;plot(SampleTs,NormalizeByRows(squeeze(SHHConvd(:,50,:)))')
end
%
tic
[MIdxs CXs Sims] = FindKepBATgAIF(DataToFit,SAIF,SHConvd,Keps1I);
tFindPK=toc;
disp(['Finding Kep, BAT took ' num2str(tFindPK)]);
TShift=MIdxs(:,1)-ZeroTIdx;
%%
% nTDif
% figure(12);clf;
% StartI=0;
% for i=1:16
%     gsubplot(16,i);
%     c=i+StartI;
%     Regressors=[HSAIF(MIdxs(c,1),:); squeeze(HSHConvd(MIdxs(c,1),Keps1I(MIdxs(c,2)),:))'];
%     CurSim=((Regressors')*CXs(:,c));
%     plot(SampleTs,DataToFit(c,:),'k',HSampleTs,CurSim,'b.',SampleTs,CurSim(1:HInterpolationFactor:end),'b*',SampleTs,DataToFit(c,:),'k*')
%     title(num2str([TShift(c) MIdxs(c,2) CXs(1,c) CXs(2,c)]));
% end
%%
PKOut=[Keps(Keps1I(MIdxs(:,2))) reshape(TDif(MIdxs(:,1)),[nToFit 1]) CXs(2,:)' CXs(1,:)'];



%% Full 3D calculation
% N=size(CTC2D,1);
% AMIdxs=zeros(N,2);
% ACXs=zeros(2,N);
% NAtATime=5000;
% 
% for i=1:NAtATime:N
%     tic
%     CurIs=i:min(N,i+NAtATime-1);
%     [AMIdxs(CurIs,:) ACXs(:,CurIs)] = FindKepBATgAIF(CTC2D(CurIs,:),SAIF,SHConvd,Keps1I);
%     t=toc;
%     disp(['Calculating ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t)]);
% end
% %%
% MeanVol=loadniidata(MeanFN);
% Kep3D=Msk2*0;
% Kep3D(Msk2)=Keps(AMIdxs(:,2));
% BAT3D=Msk2*0;
% BAT3D(Msk2)=TDif(AMIdxs(:,1));
% Vp3D=Msk2*0;
% Vp3D(Msk2)=ACXs(1,:);
% Ktrans3D=Msk2*0;
% Ktrans3D(Msk2)=ACXs(2,:);
% 
% CurSli=6;
% Tmp=squeeze(Ktrans3D(:,:,CurSli));
% Tmp2=squeeze(BAT3D(:,:,CurSli));
% Tmp2(Tmp<0.3)=0;
% figure;imagesc(mritransform(Tmp2));