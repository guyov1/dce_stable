function [PKOut OutAIFParam]=AIFTryf(WorkingP,DataToFit,MaxAmp,DataNoise,ParamAIFCoeff,nSVols,TimeBetweenDCEVols,AOptions,ShowFig,AIFFinderFN)
% DataToFit=DataToFit/1000;
TimeBetweenDCEVolsMin=TimeBetweenDCEVols/60;
InterpolationFactor=ceil(TimeBetweenDCEVols);
SampleTs=((1:nSVols)-1)*TimeBetweenDCEVolsMin;

Defaults=struct('nIterations',AOptions.EM_Num_Of_Iterations,'MaxTDif',AOptions.MaxTDif_ForAIFSearch,'SubSecRes',AOptions.SubSecondResolution);
options = optimset('TolFun',AOptions.FMS_TolFun,'MaxFunEvals',AOptions.FMS_MaxFunEvals,'MaxIter',AOptions.FMS_MaxIter,'Display','off');

Options=ExtendStruct(AOptions,Defaults);
%

nToFit=size(DataToFit,1);

T1=0.17046;A1=0.809;sig1=0.0563;sig2=0.132;A2=0.330;T2=0.365;alpha=1.050;beta=0.1685;s=38.078;tau=0.483;
tauDelta=tau-T1;T2Delta=T2-T1;
% AIFParamsA=[T1 A1 sig1];
T1=1;

% AIF_Parker10t=@(x,t) AIF_Parker( t,A1,x(3)*sig1,x(1)*T1,A2*x(5),x(3)*sig2*x(6),x(1)*T1+T2Delta*x(4),alpha*x(7),beta*x(8),s*x(9),x(1)*T1+tauDelta*x(10) )*x(2)/1000;
% AIF_Parker8t=@(x,t) AIF_Parker(
% t,A1,x(3)*sig1,x(1)*T1,A2*x(5),x(3)*sig2*x(6),x(1)*T1+T2Delta*x(4),alpha*x(7),beta*x(8),s*1,x(1)*T1+tauDelta*0 )*x(2)/1000;
% AIF_Parker8t=@(x,t) AIF_Parkerg2( t,A1,x(3)*sig1,x(1),A2*x(5),x(3)*sig2*x(6),x(1)+T2Delta*x(4),x(7),x(8))*x(2)/1000;
AIF_Parker8t=@(x,t) AIF_Parkerg2( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8))*x(2);
MinFirstBolusSig=AOptions.MinFirstBolusStd; %2; % seconds
LB=[0.1 0     MinFirstBolusSig/60 0.1     0   0.1  0.1 0]';
UB=[10  1.5   0.25                10      1   3    2   0.3]';
% X0=[T1  A1  sig1                T2Delta A2  sig2 Alpha Beta]';
X0=[1   1     0.15                1 0.3 0.3   0.4 0.2]';
HInterpolationFactor=ceil(InterpolationFactor*Options.SubSecRes);
Hdt=TimeBetweenDCEVolsMin/HInterpolationFactor;
HSampleTs=0:Hdt:SampleTs(end);
nHVols=numel(HSampleTs);
% nKeps=100;
% Keps=gpowspace(0,15,nKeps,5)';
% nKeps1=100;
% Keps1I=floor(linspace(1,nKeps,nKeps1));
WReg=3;
WRegBolus=3;
WSmooth1=30;
WSmooth1Bolus=10;
WSmooth2=0.000;
WSmooth2Bolus=0.000;
% TDif=(-0*TimeBetweenDCEVolsMin):Hdt:(2*TimeBetweenDCEVolsMin);
ThreeSec=ceil(Options.MaxTDif/(Hdt*60));
TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;
nTDif=numel(TDif);
ZeroTIdx=find(TDif==0);
BaseParams=ParamAIFCoeff(1:8,3);
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
% Start with AIF from Artery
% VesselsIdxs=find(BinCVI==max(BinCVI));
VesselsIdxs=nToFit;

% MaxAmp=max(max(DataToFit(VesselsIdxs,:)));
% MaxAmp=max(max(DataToFit));
NormalizedVessels=NormalizeByRows(DataToFit(VesselsIdxs,:))*MaxAmp;
MeanVessel=mean(NormalizedVessels,1);
AIF_Parker8tx=@(x,t) AIF_Parker8t(x,t).*MaxAmp;
[OldParams1,resnorm,residual,exitflag,output] = lsqcurvefit(AIF_Parker8tx,X0',SampleTs,MeanVessel,LB, UB,options);

HAIF=AIF_Parker8t(OldParams1,HSampleTs);
OldParams=OldParams1';
OldParams(2)=OldParams(2)*MaxAmp;
OldParamsFromMeanVessel=OldParams;
UB(2)=UB(2)*MaxAmp;
HAIF=AIF_Parker8t(OldParams,HSampleTs);
HAIFFromMeanVessel=AIF_Parker8t(OldParamsFromMeanVessel,HSampleTs);

BAIF=AIF_Parker8tx(OldParams1,HSampleTs);
X0AIF=AIF_Parker8tx(X0,HSampleTs);
% figure(99);plot(SampleTs,MeanVessel,'k-*');hold on;plot(HSampleTs,HAIF,'b','LineWidth',3);
figure(1000);clf;
if(Options.nIterations~=0)
    subplot(2,2,1);
end
hold on;plot(HSampleTs,X0AIF,'r','LineWidth',2);plot(HSampleTs,HAIFFromMeanVessel,'b','LineWidth',2);plot(HSampleTs,BAIF,'g','LineWidth',1);plot(SampleTs,MeanVessel,'k*');
% figure(98);plot(SampleTs,AIF_Parker8tx(OldParams1,SampleTs),'g',HSampleTs,AIF_Parker8tx(OldParams1,HSampleTs),'b')

% X0=[T1  A1  sig1 T2Delta A2 sig2 Alpha Beta]';
% X0=[1  1  0.15 1 0.3 0.3   0.4 0.2]';
% figure;plot(SampleTs,MeanVessel,'k-*',SampleTs,AIF_Parker8tx(X0,SampleTs),'g')
%%
WPerCTC=1./(max(DataToFit,[],2).*(DataNoise.^0));
WPerCTC=10000*nHVols*WPerCTC./(sum(WPerCTC).*nSVols);    
    
HBolusStart=ceil(OldParams(1)/Hdt);
nTimePoints=numel(SampleTs);
nHTimePoints=numel(HSampleTs);
ConvIdxM=CreateConvIdxMFromSampleTs(numel(SampleTs));
TriB=ConvIdxM>0;
% ConvIdxMTriB=ConvIdxM(TriB);
HConvIdxM=CreateConvIdxMFromSampleTs(numel(HSampleTs));
HTriB=HConvIdxM>0;
HConvIdxMTriB=HConvIdxM(HTriB);

if(Options.nIterations==0)
    OutAIFParamBeforeFMS=OldParams';
else
    tic
    % HBolusStart=(BolusStart-1)*HInterpolationFactor+1;

    HAIF=AIF_Parker8t(OldParams,HSampleTs);
    tBase=toc;
    disp(['Preparation took ' num2str(tBase)]);
    %%
    IterCont=true;
    IterCount=0;
    clear AllAIF2s AllAIF2cs AllAIFPs AllRMSs ARMS
    AllAIF2cs(1,:)=HAIF;
    BeforeIterAIF=HAIF;
    AllAIFPs(1,:)=OldParams;
    while(IterCont)
        IterCount=IterCount+1;
        IterCont=IterCount<Options.nIterations;
        disp(['Iteration #' num2str(IterCount)]);
        %
        tic;
        %     HConvd=DCECostFuncgrT1ForConv(HAIF',Keps,HSampleTs,HConvIdxMTriB,HTriB);
        tConv=toc;
        %     disp(['Convolving the AIF took ' num2str(tConv)]);
        %     if(ShowFig)
        %         figure;plot(HSampleTs,NormalizeByRows(HConvd(:,:))')
        %     end
        %
        tic
        CHAIF=cumtrapz(HSampleTs,HAIF);
        SAIF=zeros([nTDif numel(SampleTs)]);
        CSAIF=zeros([nTDif numel(SampleTs)]);
        for i=1:nTDif
            SAIF(i,:)=interp1(HSampleTs,HAIF,SampleTs+TDif(i),[],'extrap');
            CSAIF(i,:)=interp1(HSampleTs,CHAIF,SampleTs+TDif(i),[],'extrap');
        end
        tSampleConv=toc;
        %     disp(['Sampling convolved CTCs took ' num2str(tSampleConv)]);
        %     if(ShowFig)
        %         figure;plot(SampleTs,NormalizeByRows(squeeze(SHHConvd(:,50,:)))')
        %     end
        %
        tic
        %     [MIdxs CXs Sims] = FindKepBATgAIF(DataToFit,SAIF,SHConvd,Keps1I);
        CumSumAIFdTSamp=CSAIF;
        [PKs Sims ARMS(IterCount)] = FindPKBATgAIFMurase(WorkingP,DataToFit,SAIF,SampleTs,CumSumAIFdTSamp,WPerCTC);
        %     CXs=PKs(:,[3 4])';
        %     MIdxs=PKs(:,1);

        tFindPK=toc;
        %     disp(['Finding Kep, BAT took ' num2str(tFindPK)]);
        %     ARMS(IterCount)=sqrt(mean(((Sims-DataToFit).^2),2))'*WPerCTC;
        disp(['Arms - ' num2str(ARMS(IterCount))]);
        % Deconv
        tic
        RHS=[];
        LHS=[];

        % BAT Kep Vp Ktrans Ve RMS
        [U, QQQ, IB]=unique(PKs(:,2));
        nUKep=numel(U);
        EKepM=zeros(nHTimePoints);
        EKepMs=zeros([nUKep nTimePoints nHTimePoints]);
        for i=1:nUKep
            EKepV=Hdt*exp(-(U(i))*HSampleTs);
            EKepM(HTriB)=EKepV(HConvIdxMTriB);
            EKepMs(i,:,:)=EKepM(1:HInterpolationFactor:end,:); %
        end
        TShift=PKs(:,1)-ZeroTIdx;
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

            CurLHS=EKepM*PKs(c,4)+LittleI*PKs(c,3);
            GoodRows=sum(LittleI,2)>0;
            CurLHS=CurLHS*WPerCTC(c);

            LHS=[LHS; CurLHS(GoodRows,:)];
            CurRHS=DataToFit(c,:)'*WPerCTC(c);

            RHS=[RHS; CurRHS(GoodRows,:)];
            %         figure;plot(CurRHS(GoodRows,:));hold on;plot(CurLHS(GoodRows,:)*(HAIF'),'g')
        end
        % Regularization
        RefAIF=AllAIF2cs(max(IterCount-1,1),:)';
        %     RefAIF=HAIF';
        AroundBolus=(1:nHTimePoints);
        AroundBolus=AroundBolus<HBolusStart+HInterpolationFactor*1 & AroundBolus>HBolusStart-HInterpolationFactor*1;
        LHSReg=diag(AroundBolus*WRegBolus+(~AroundBolus)*WReg);
        RHSReg=RefAIF*WReg;
        LHSSmooth1=diag(AroundBolus*WSmooth1Bolus+(~AroundBolus)*WSmooth1);
        LHSSmooth1=LHSSmooth1(1:end-1,:)-LHSSmooth1(2:end,:);
        RHSSmooth1=zeros(nHTimePoints-1,1)+1*LHSSmooth1*RefAIF;%  zeros(nHTimePoints-1,1);
        if(WSmooth2+WSmooth2Bolus>0)
            LHSSmooth2=diag(AroundBolus*WSmooth2Bolus+(~AroundBolus)*WSmooth2);
            LHSSmooth2=LHSSmooth2(1:end-2,:)-2*LHSSmooth2(2:end-1,:)+LHSSmooth2(3:end,:);
            RHSSmooth2=LHSSmooth2*RefAIF;
            LHS=[LHS; LHSReg; LHSSmooth1; LHSSmooth2];
            RHS=[RHS; RHSReg; RHSSmooth1; RHSSmooth2];
        else
            LHS=[LHS; LHSReg; LHSSmooth1];
            RHS=[RHS; RHSReg; RHSSmooth1];
        end
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
        %     NewPc = lsqcurvefit(AIF_Parker8t,OldParams(1:8),HSampleTs,AIF2',ParamAIFCoeff(1:8,1), ParamAIFCoeff(1:8,2),options)';

        [NewPc,resnorm,residual,exitflag,outputs(IterCount)] = lsqcurvefit(AIF_Parker8t,OldParams(1:8),HSampleTs,AIF2',LB, UB,options);

        %     x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options)
        NewPc(2)=OldParams(2);
        NewPc(1)=OldParams(1);
        AIF2c=AIF_Func1(NewPc)';
        RMS2=gCost(LHS*AIF2c,RHS,'RMS');
        tToFuncParams=toc;
        disp(['Converting to functional form took ' num2str(tToFuncParams) ' to RMS ' num2str(RMS2)]);
        %     if(ShowFig)
        %         figure;plot(HSampleTs,[HAIF' AIF2 AIF2c])
        %     end
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
    %
    [Tmp, MinARMSI]=min(ARMS)
    MinARMSI=max(MinARMSI,2);
    % MinARMSI=1;
    %
    % if(ShowFig)
    subplot(2,2,2);
    plot(AllAIF2s');

    subplot(2,2,3);
    plot(HSampleTs,AllAIF2cs');hold on;
    plot(SampleTs,MeanVessel,'r-*','LineWidth',4);
    plot(HSampleTs,BeforeIterAIF,'k','LineWidth',2);
    plot(HSampleTs,AllAIF2cs(MinARMSI-1,:),'b','LineWidth',2);
    subplot(2,2,4);
    plot([AllRMSs ARMS'/nToFit]);
    title(MinARMSI);
    % end
    OutAIFParamBeforeFMS=AllAIFPs(MinARMSI,:);
end
%%
optionsB=options;
optionsB.Display='final';
optionsB.TolFun=10^-3;
RestrainX=@(x) min(UB',max(LB',x));
% CostFunc6=@(x) -0*min(x(1),UB(3))/8 -0*abs(min(UB(5),max(LB(5),x(3)))-0.3)/20 +AIFCostFuncM(WorkingP,RestrainX([AllAIFPs(1,1:2) x(1:6)]),DataToFit,SampleTs,HSampleTs,TDif,WPerCTC);
CostFunc6=@(x) -0*min(x(1),UB(3))/8 -0*abs(min(UB(5),max(LB(5),x(3)))-0.3)/20 +AIFCostFuncM(WorkingP,RestrainX([OldParams(1:2)' x(1:6)]),DataToFit,SampleTs,HSampleTs,TDif,WPerCTC);
% BaseParams=AllAIFPs(2,:);
BaseParams=OutAIFParamBeforeFMS;
% BaseParams=X0';
% X0=[T1  A1  sig1 T2Delta A2 sig2 Alpha Beta]';
x = fminsearch(CostFunc6,BaseParams(3:8),optionsB)
FMSCost=CostFunc6(x);
% OutAIFParam=RestrainX([AllAIFPs(1,1:2) x(1:6)]);
OutAIFParam=RestrainX([OldParams(1:2)' x(1:6)]);
FMSAIF=AIF_Parker8t(OutAIFParam,HSampleTs);
if(Options.nIterations~=0)
    figure(1000);
    subplot(2,2,3);hold on;
end
plot(HSampleTs,FMSAIF,'m','LineWidth',2);
% IterAIF=AIF_Parker8t(OutAIFParamBeforeFMS,HSampleTs);
% plot(IterAIF,'g','LineWidth',2);
[UB';LB';OutAIFParamBeforeFMS;OutAIFParam]
%%
gprint(1000,[AIFFinderFN(1:end-4) '.jpg']);
close(1000);
%%
HAIF=AIF_Parker8t(OutAIFParam,HSampleTs);
CHAIF=cumtrapz(HSampleTs,HAIF);
CHAIF=[0 cumsum(HAIF(1:end-1))]*Hdt;

SAIF=zeros([nTDif numel(SampleTs)]);
CSAIF=zeros([nTDif numel(SampleTs)]);
for i=1:nTDif
    SAIF(i,:)=interp1(HSampleTs,HAIF,SampleTs+TDif(i),[],'extrap');
    CSAIF(i,:)=interp1(HSampleTs,CHAIF,SampleTs+TDif(i),[],'extrap');
end
[PKs Sims] = FindPKBATgAIFMuraseF(DataToFit,SAIF,SampleTs,CSAIF);

% [PKs Sims] = FindPKBATgAIFMurase(WorkingP,DataToFit,SAIF,SampleTs,CSAIF);
CostsA=sqrt(mean(((Sims-DataToFit).^2),2));
Costs=CostsA.*WPerCTC;
%%
nKeps=100;
Keps=gpowspace(0,15,nKeps,5)';
nKeps1=100;
Keps1I=floor(linspace(1,nKeps,nKeps1));

HConvd=DCECostFuncgrT1ForConv(HAIF',Keps,HSampleTs,HConvIdxMTriB,HTriB);
%
SHConvd=zeros([nTDif nKeps numel(SampleTs)]);
HSHConvd=zeros([nTDif nKeps numel(HSampleTs)]);
for i=1:nTDif
    SHConvd(i,:,:)=interp1(HSampleTs,HConvd',SampleTs+TDif(i),[],'extrap')';
    HSHConvd(i,:,:)=interp1(HSampleTs,HConvd',HSampleTs+TDif(i),[],'extrap')';
end
%
tic
[MIdxs CXs Sims] = FindKepBATgAIF(DataToFit,SAIF,SHConvd,Keps1I);
tFindPK=toc;
disp(['Finding Kep, BAT took ' num2str(tFindPK)]);
%%
CurKeps=PKs(:,2);
% CurKeps(18)=0;
HConvd2=DCECostFuncgrT1ForConv(HAIF',CurKeps,HSampleTs,HConvIdxMTriB,HTriB);
%
tic
HSHConvd2=zeros([nTDif nToFit numel(HSampleTs)]);
HSAIF=zeros([nTDif numel(HSampleTs)]);
for i=1:nTDif
    SHConvd2(i,:,:)=interp1(HSampleTs,HConvd2',SampleTs+TDif(i),[],'extrap')';
    HSHConvd2(i,:,:)=interp1(HSampleTs,HConvd2',HSampleTs+TDif(i),[],'extrap')';
    HSAIF(i,:)=interp1(HSampleTs,HAIF,HSampleTs+TDif(i),[],'extrap');
end
tSampleConv=toc;
disp(['Sampling convolved CTCs took ' num2str(tSampleConv)]);
% if(ShowFig)
%     figure;plot(SampleTs,NormalizeByRows(squeeze(SHHConvd(:,50,:)))')
% end
%
TShift=PKs(:,1)-ZeroTIdx;
% if(ShowFig)
    figure(12);clf;
    StartI=0;
    for i=1:nToFit
        gsubplot(nToFit,i);
        c=i+StartI;
        Regressors=[SAIF(PKs(c,1),:); squeeze(SHConvd2(PKs(c,1),c,:))'];
        % BAT Kep Vp Ktrans Ve RMS fVal
        CurSim2=((Regressors')*(PKs(c,[3 4])'));
        Regressors=[HSAIF(MIdxs(c,1),:); squeeze(HSHConvd(MIdxs(c,1),Keps1I(MIdxs(c,2)),:))'];
        CurSim=((Regressors')*(CXs(:,c)));
        plot(SampleTs,DataToFit(c,:),'k',HSampleTs,CurSim,'b.',SampleTs,CurSim(1:HInterpolationFactor:end),'b*',SampleTs,DataToFit(c,:),'k*')
        hold on; plot(SampleTs,Sims(c,:),'r');plot(SampleTs,CurSim2,'g');
        RError=sqrt(mean((CurSim2'-DataToFit(c,:)).^2))/max(DataToFit(c,:));
        MError=sqrt(PKs(c,7)/nToFit)/max(DataToFit(c,:));
        %             title(['  ' num2str([TShift(c) PKs(c,2) PKs(c,3) PKs(c,4)])]);
        %     title(['  ' num2str([c CostsA(c) Costs(c) WPerCTC(c)])]);
        %         title(PKs(c,7));
        %         title([RError MError]);
        title([max(CurSim)./max(DataToFit(c,:)) TShift(c)]) ;
        set(gca,'YTick',[]);
        set(gca,'XTick',[]);
        HSims(i,:)=CurSim;
    end
    saveas(12,[AIFFinderFN(1:end-4) '_FitOnArt.fig']);
    gprint(12,[AIFFinderFN(1:end-4) '_FitOnArt.jpg']);
    close(12);
% end
%%
if(nargin>8)
    if(Options.nIterations==0)
        save(AIFFinderFN,'OutAIFParam');
    else
        save(AIFFinderFN,'AllAIF2s','AllAIF2cs','AllRMSs','ARMS','AllAIFPs','OutAIFParam');
    end
end
% PKOut=[Keps(Keps1I(MIdxs(:,2))) reshape(TDif(MIdxs(:,1)),[nToFit 1]) CXs(2,:)' CXs(1,:)'];
PKOut=[PKs(:,2) reshape(TDif(PKs(:,1)),[nToFit 1]) PKs(:,4) PKs(:,3)];