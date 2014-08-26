function [PKOut OutAIFParam HAIF]=AIFTryf(WorkingP,DataToFit,MaxAmp,DataNoise,ParamAIFCoeff,nSVols,TimeBetweenDCEVols,AOptions,ShowFig,AIFFinderFN,TimeVec,Handle)
% DataToFit=DataToFit/1000;

%% Initial params

% Get the between volumes and create time samples vector
TimeBetweenDCEVolsMin=double(TimeBetweenDCEVols)/60;
InterpolationFactor=ceil(TimeBetweenDCEVols);
SampleTs=((1:nSVols)-1)*TimeBetweenDCEVolsMin;

nSVols=numel(TimeVec);
SampleTs=TimeVec;
DTimeVec=diff(TimeVec);
F=find(DTimeVec>DTimeVec(1)*2,1);
nNormalTs=numel(TimeVec);
if(~isempty(F))
    nNormalTs=F(1);
end

Defaults=struct('nIterations',AOptions.EM_Num_Of_Iterations,'MaxTDif',AOptions.MaxTDif_ForAIFSearch,'SubSecRes',AOptions.SubSecondResolution);
options = optimset('TolFun',AOptions.FMS_TolFun,'MaxFunEvals',AOptions.FMS_MaxFunEvals,'MaxIter',AOptions.FMS_MaxIter,'Display','off');

% Extend Options (and not options!) struct with the Defaults just determined
Options=ExtendStruct(AOptions,Defaults);
%

%% Setting the super resolution parameters

% Number of representing voxels to fit
nToFit=size(DataToFit,1);
% Parameters of representing AIF (like Parker) I assume (should calrify with Gilad)
T1=0.17046;A1=0.809;sig1=0.0563;sig2=0.132;A2=0.330;T2=0.365;alpha=1.050;beta=0.1685;s=38.078;tau=0.483;
tauDelta=tau-T1;T2Delta=T2-T1;
% AIFParamsA=[T1 A1 sig1];
T1=1;

% AIF_Parker10t=@(x,t) AIF_Parker( t,A1,x(3)*sig1,x(1)*T1,A2*x(5),x(3)*sig2*x(6),x(1)*T1+T2Delta*x(4),alpha*x(7),beta*x(8),s*x(9),x(1)*T1+tauDelta*x(10) )*x(2)/1000;
% AIF_Parker8t=@(x,t) AIF_Parker(
% t,A1,x(3)*sig1,x(1)*T1,A2*x(5),x(3)*sig2*x(6),x(1)*T1+T2Delta*x(4),alpha*x(7),beta*x(8),s*1,x(1)*T1+tauDelta*0 )*x(2)/1000;
% AIF_Parker8t=@(x,t) AIF_Parkerg2( t,A1,x(3)*sig1,x(1),A2*x(5),x(3)*sig2*x(6),x(1)+T2Delta*x(4),x(7),x(8))*x(2)/1000;

% Create a function handle to AIF_Parkerg2 which calculates the AIF according to Parker
% C=AIF_Parkerg2(t,A1,sig1,T1,rA2,sig2,T2,ralpha,beta)
AIF_Parker8t=@(x,t) AIF_Parkerg2( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8))*x(2);
AIF_Parker9t=@(x,t) AIF_Parkerg3( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8),x(9),max(TimeVec))*x(2);

% Taken from user options in GUI
%The minimum width of the bolus (standard deviation of the Gaussian that represents the first bolus).
MinFirstBolusSig=AOptions.MinFirstBolusStd; %2; % seconds

% ASK GILAD - Explain the meaning of all the following numbers...
% ASWER - 
LB=[0.1 0     MinFirstBolusSig/60 0.1     0   0.1  0.0001   0  0]';
UB=[10  1.5   0.25                10      1   1     2        2  1]';

% X0=[T1  A1  sig1                T2Delta A2  sig2 Alpha Beta BaseLevel]';
X0=[1   1     0.15                1       0.3 0.3   0.4  0.2  0.1]';

% Setting the super resolution to be every 0.5 seconds (or ~ 0.0083 minute)
HInterpolationFactor=ceil(InterpolationFactor*Options.SubSecRes);
Hdt=TimeBetweenDCEVolsMin/HInterpolationFactor;
% HSampleTs=0:Hdt:SampleTs(end);
HSampleTs=[0:Hdt:SampleTs(nNormalTs) TimeVec((nNormalTs+1):end)];
nHVols=numel(HSampleTs);

% Possible Keps values
% nKeps=100;

% ASK GILAD - What is the meaning of the gpow space instead of the linear one? that the intervals are not identical?
% ANSWER - He wanted a bigger resolution in the smaller numbers. Linear space would not necessarily help when trying to find the best kep.
% Get 100 possible values for Kep in the range 0-15
% Keps=gpowspace(0,15,nKeps,5)';
% nKeps1=100;
% Set linear space from 1-100 (1,2,3,...100)  - > Index of pssible Keps
% Keps1I=floor(linspace(1,nKeps,nKeps1));

% Some more parameters
WReg=3;
WRegBolus=3;
WSmooth1=30;
WSmooth1Bolus=10;
WSmooth2=0.000;
WSmooth2Bolus=0.000;
% TDif=(-0*TimeBetweenDCEVolsMin):Hdt:(2*TimeBetweenDCEVolsMin);

% Number of super resolution time stams in 3 seconds
ThreeSec=ceil(Options.MaxTDif/(Hdt*60));

% Create a linear space about -3:0.5:+3 seconds
TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;
% Number of time stamps in the -3 +3 interval
nTDif=numel(TDif);

% Zero time index
ZeroTIdx=find(TDif==0);

% Get the parameters calculated before in DCET1_PKf
%OldParams=ParamAIFCoeff(1:8,3);
BaseParams=ParamAIFCoeff(1:8,3);

%
%  (lower bound, upper bound and start point)
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

% Change some of the parameters bounds
ParamAIFCoeff(3,1)=1;
ParamAIFCoeff(2,[1 2])=[0.01 1000];


%% Start with AIF from  Artery

% Get the voxels from the maximum cluster -> highest
% ASK GILAD - Why does he decide these are artery voxels? Don't all the voxels in the tumor get a peak, same as the ones in artery?
% ANSWER - BinCVI are currently the voxels with high value around bolus. It will change to ones eith also low value around the end.
% VesselsIdxs=find(BinCVI==max(BinCVI));

VesselsIdxs=nToFit;

% Maximum Amplitude
% MaxAmp=max(max(DataToFit(VesselsIdxs,:)));
% MaxAmp=max(max(DataToFit));

% NormalizeByRows first divides each c(t) by its maximal value.
% We then multiply by the general maximal value.
% This way all C(t)s are normalized according to the maximal amplitude of all C(t)s
NormalizedVessels=NormalizeByRows(DataToFit(VesselsIdxs,:))*MaxAmp;

% Get the mean c(t) of the vessel
MeanVessel=mean(NormalizedVessels,1);

% Create a function handle which normalizes Parker's AIF according to max amplitude
AIF_Parker8tx=@(x,t) AIF_Parker8t(x,t).*MaxAmp;
AIF_Parker9tx=@(x,t) AIF_Parker9t([x(1:8) x(9)/MaxAmp],t).*MaxAmp;
% ASK GILAD - What parameters does he provide the least square fit curve function? and how does it work exactly?
% ANSWER - lsqcurvefit will minimize the squared error and will return the
% AIF_Parker8tx parameters that minimize it.
% He gives it initial parameters ( OldParams(1:8) ), time domain, mean voxel c(t) and lower and upper bound for the 8 paramters.
% [OldParams1,resnorm,residual,exitflag,output] = lsqcurvefit(AIF_Parker8tx,X0',SampleTs,MeanVessel,LB, UB,options);

% CostFuncMeanVessel=@(x) gCost(MeanVessel,AIF_Parker9tx(min(max(x,[LB; 0]'),[UB; 1]'),SampleTs),'RMS');
CostFuncMeanVessel=@(x) gCost(MeanVessel,AIF_Parker9tx(min(max(x,LB'),UB'),SampleTs),'RMS');
% OldParams1=fminsearch(CostFuncMeanVessel,[X0(1:8)' MeanVessel(end)]);
X0(end)=MeanVessel(end);
OldParams1=fminsearch(CostFuncMeanVessel,[X0']);
OldParams1=min(max(OldParams1,LB'),UB');

% figure(121);clf;
% plot(SampleTs,[MeanVessel; AIF_Parker9tx(OldParams1,SampleTs)],'.-')

% Parker's AIF using the parameters calculated by the curve fitting
HAIF=AIF_Parker9tx(OldParams1,HSampleTs);
% gfig;plot(SampleTs,MeanVessel,'k*');hold on;plot(HSampleTs,HAIF);
% hold on;
% plot(HSampleTs,HAIF,'r.-')
% ASK GILAD - Why did he use the following normalization? What is it good for?
% ANSWER - He wanted to normalize it according to the maximal value he gotfor all c(t). 
%                    This is why he divided it by the maximal value it just had and multiplied by the new value he wants.
%                    Multiply  by the orignal Max Amplitude and divide by the new one created by Parker's method. Then, recalculate the AIF.
OldParams=OldParams1';
OldParams(2)=OldParams(2)*MaxAmp;
OldParams(9)=OldParams(9)/MaxAmp;
OldParamsFromMeanVessel=OldParams;
UB(2)=UB(2)*MaxAmp;

% Calculate Parker's AIF again using the new parameters
HAIF=AIF_Parker9tx(OldParams',HSampleTs);
HAIFFromMeanVessel=AIF_Parker9t(OldParamsFromMeanVessel',HSampleTs);

BAIF=AIF_Parker9tx(OldParams1,HSampleTs);
X0AIF=AIF_Parker8tx(X0,HSampleTs);
% figure(99);plot(SampleTs,MeanVessel,'k-*');hold on;plot(HSampleTs,HAIF,'b','LineWidth',3);

figure(Handle);clf;
if(Options.nIterations~=0)
    subplot(2,2,1);
end

hold on;plot(HSampleTs,X0AIF,'r','LineWidth',2);plot(HSampleTs,HAIFFromMeanVessel,'b','LineWidth',2);plot(HSampleTs,BAIF,'g','LineWidth',1);plot(SampleTs,MeanVessel,'k*');
title('r-Start AIF. b,k,g - Mean vessel AIF.');
% figure(98);plot(SampleTs,AIF_Parker8tx(OldParams1,SampleTs),'g',HSampleTs,AIF_Parker8tx(OldParams1,HSampleTs),'b')

% X0=[T1  A1  sig1 T2Delta A2 sig2 Alpha Beta]';
% X0=[1  1  0.15 1 0.3 0.3   0.4 0.2]';
% figure;plot(SampleTs,MeanVessel,'k-*',SampleTs,AIF_Parker8tx(X0,SampleTs),'g')

%%

% Taking each representing voxel's maximal value and relative noise and calculate for each:
% 1 ./ ( max_val .* r_noise ) ->  1 / noise which will be the weight of that voxel
WPerCTC=1./(max(DataToFit,[],2).*(DataNoise.^0));

WPerCTC=10000*nHVols*WPerCTC./(sum(WPerCTC).*nSVols);    
    
% Get the Bolus arrival time stamp index (volume number in super resolution)
HBolusStart=ceil(OldParams(1)/Hdt);

% Original number of time stamps
nTimePoints=numel(SampleTs);

% Super resolution number of time stamps
nHTimePoints=numel(HSampleTs);

% Create convolution indices. Example when there are 50 time stamps:
% 1   0   0   ... 0
% 2   1   0   ... 0
% 3   2   1   ... 0
%      ...
% 49 48 47 ... 0
% 50 49 48 ... 1
ConvIdxM=CreateConvIdxMFromSampleTs(numel(SampleTs));

% Mask for all indices bigger than 0
TriB=ConvIdxM>0;

% ConvIdxMTriB=ConvIdxM(TriB);

% Do the same as above for the super resolution
HConvIdxM=CreateConvIdxMFromSampleTs(numel(HSampleTs));

% Mask for all indices bigger than 0
HTriB=HConvIdxM>0;

% Get all indices with value > 0 column vice (   this is why we will see in the beginning [1 ; 2 ; 3 ... ]   )
HConvIdxMTriB=HConvIdxM(HTriB);

if(nNormalTs~=numel(SampleTs))
    Options.nIterations=0;
end
% If user chose 0 iterations, don't calculate the EM part (only the FMS)
if(Options.nIterations==0)
    OutAIFParamBeforeFMS=OldParams';
else
    %% Parker's AIF
    tic
    % HBolusStart=(BolusStart-1)*HInterpolationFactor+1;
    
    % Calculate Parker's AIF again using the new parameters
    HAIF=AIF_Parker8t(OldParams,HSampleTs);
    
    tBase=toc;
    disp(['Preparation (calculation of Parkers AIF took ' num2str(tBase)]);
    
    %% Initiate iterations parameters
    IterCont=true;
    IterCount=0;
    clear AllAIF2s AllAIF2cs AllAIFPs AllRMSs ARMS
    
    % Concatenation of all AIFs
    AllAIF2cs(1,:)=HAIF;
    % AIF before the iteration
    BeforeIterAIF=HAIF;
    % Concatenation of all AIF params
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
        
        % Calculate the commulutive intergral of the AIF
        CHAIF=cumtrapz(HSampleTs,HAIF);
        
        SAIF=zeros([nTDif numel(SampleTs)]);
        CSAIF=zeros([nTDif numel(SampleTs)]);
        
        
        for i=1:nTDif
            
            % Interpolate AIF(t) over the original time points shifted by TDif(i)
            SAIF(i,:)=interp1(HSampleTs,HAIF,SampleTs+TDif(i),[],'extrap');
            
            % Interpolate the cummulative intergral of the AIF(t) over the original time points shifted by TDif(i)
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
        AllAIFPs(IterCount+1,1:8)=NewPc;
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
    title('AIFs from EM iterations');

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
    title('EM iteration scores');
end
%%
optionsB=options;
optionsB.Display='final';
optionsB.Display='iter';
optionsB.TolFun=10^-3;
% RestrainX=@(x) min(UB',max(LB',x));
% CostFunc6=@(x) -0*min(x(1),UB(3))/8 -0*abs(min(UB(5),max(LB(5),x(3)))-0.3)/20 +AIFCostFuncM(WorkingP,RestrainX([AllAIFPs(1,1:2) x(1:6)]),DataToFit,SampleTs,HSampleTs,TDif,WPerCTC);
% CostFunc6=@(x) -0*min(x(1),UB(3))/8 -0*abs(min(UB(5),max(LB(5),x(3)))-0.3)/20 +AIFCostFuncM(WorkingP,RestrainX([OldParams(1:2)' x(1:6)]),DataToFit,SampleTs,HSampleTs,TDif,WPerCTC);

% X0=[T1  A1  sig1                T2Delta A2  sig2 Alpha Beta BaseLevel]';
% UB(end)=OutAIFParamBeforeFMS(end);
% LB(end)=OutAIFParamBeforeFMS(end);

Idx4=[1 5 6 7];
% RestrainX=@(x) min([UB; 1]',max([LB; 0]',x));
RestrainX3=@(x) min(UB(3:end)',max(LB(3:end)',x));
RestrainX=@(x) min(UB',max(LB',x));
RestrainX4=@(x) min(UB(Idx4+2)',max(LB(Idx4+2)',x));
% SimilarityCost3=@(x) gCost(OutAIFParamBeforeFMS(3:9),RestrainX3(x),'SumAbs');
SimilarityCost3=@(x) gCost(OutAIFParamBeforeFMS(Idx4+2),RestrainX4(x(Idx4)),'SumAbs');
ShapeSimCost7=@(x) abs(AIF_Parker9t(OutAIFParamBeforeFMS,SampleTs(end))-AIF_Parker9t(RestrainX([OldParams(1:2)' x(1:7)]),SampleTs(end)));
ShapeSimCostb7=@(x) gCost(AIF_Parker9t(OutAIFParamBeforeFMS,SampleTs),AIF_Parker9t(RestrainX([OldParams(1:2)' x(1:7)]),SampleTs),'RMS');
FitCost7=@(x) AIFCostFuncMf(WorkingP,RestrainX([OldParams(1:2)' x(1:7)]),DataToFit,SampleTs,HSampleTs,TDif,WPerCTC,AIF_Parker9t);
CostFunc7=@(x) FitCost7(x)+SimilarityCost3(x)*Options.WeightForAIFMeanVesses+ShapeSimCost7(x)*1+ShapeSimCostb7(x)*1000;
% CostFunc7=@(x) -0*min(x(1),UB(3))/8 -0*abs(min(UB(5),max(LB(5),x(3)))-0.3)/20+ + +AIFCostFuncMf(WorkingP,RestrainX([OldParams(1:2)' x(1:7)]),DataToFit,SampleTs,HSampleTs,TDif,WPerCTC,AIF_Parker9t);
SimilarityCost4=@(x) gCost(OutAIFParamBeforeFMS(Idx4+2),RestrainX4(x),'SumAbs');
ShapeSimCost4=@(x) abs(AIF_Parker9t(OutAIFParamBeforeFMS,SampleTs(end))-AIF_Parker9t(RestrainX([OldParams(1:2)' x(1) 1 0 1 x(2:4)]),SampleTs(end)));
ShapeSimCostb4=@(x) gCost(AIF_Parker9t(OutAIFParamBeforeFMS,SampleTs),AIF_Parker9t(RestrainX([OldParams(1:2)' x(1) 1 0 1 x(2:4)]),SampleTs),'RMS');
FitCost4=@(x) AIFCostFuncMf(WorkingP,RestrainX([OldParams(1:2)' x(1) 1 0 1 x(2:4)]),DataToFit,SampleTs,HSampleTs,TDif,WPerCTC,AIF_Parker9t);
CostFunc4=@(x) FitCost4(x)+SimilarityCost4(x)*Options.WeightForAIFMeanVesses+ShapeSimCost4(x)*1+ShapeSimCostb4(x)*1000;

swarn=warning('off','MATLAB:interp1:EmptyMethod');
% BaseParams=AllAIFPs(2,:);
BaseParams=OutAIFParamBeforeFMS;
x = fminsearch(CostFunc4,BaseParams(Idx4+2),optionsB)
Base4Params=[OldParams(1:2)' 1 1 0 1 2 3 4];
Best4Paramsa=InsertIntoIdxs(Base4Params,Idx4+2,x);
Best4Params=RestrainX([OldParams(1:2)' Best4Paramsa(3:9)]);
% gfig(2102);plot(AIF_Parker9t(Best4Params,HSampleTs));

% BaseParams=X0';
% X0=[T1  A1  sig1 T2Delta A2 sig2 Alpha Beta]';
% x = fminsearch(CostFunc7,BaseParams(3:9),optionsB)
x = fminsearch(CostFunc7,Best4Params(3:9),optionsB)
% x = OldParamsFromMeanVessel(3:9)';
FMSCost=CostFunc7(x);

% OutAIFParam=RestrainX([AllAIFPs(1,1:2) x(1:6)]);
OutAIFParam=RestrainX([OldParams(1:2)' x(1:7)]);
Best7Params=OutAIFParam;
FMSFCostMean=FitCost7(OutAIFParamBeforeFMS(3:end));
FMSFCost7=FitCost7(OutAIFParam(3:end));
FMSFCost4=FitCost7(Best4Params(3:end));
FMSCostMean=CostFunc7(OutAIFParamBeforeFMS(3:end));
FMSCost7=CostFunc7(OutAIFParam(3:end));
FMSCost4=CostFunc7(Best4Params(3:end));
[FMSFCostMean FMSFCost7 FMSFCost4]
[FMSCostMean FMSCost7 FMSCost4]
if(FMSCost4<FMSCost7)
    OutAIFParam=Best4Params;
end
[OutAIFParamBeforeFMS; OutAIFParam; Best4Params]
FMSAIF=AIF_Parker9t(OutAIFParam,HSampleTs);
% figure;plot([AIF_Parker9t(OutAIFParamBeforeFMS,HSampleTs); AIF_Parker9t(OutAIFParam,HSampleTs)]')
if(Options.nIterations~=0)
    figure(Handle);
    subplot(2,2,3);hold on;
    plot(HSampleTs,FMSAIF,'m','LineWidth',2);
    title('r,k-mean Vessel, b-after EM, m-Chosen AIF');
else
    plot(HSampleTs,FMSAIF,'m','LineWidth',2);
    title('r-Start AIF. b,k,g - Mean vessel AIF, m-Chosen AIF');
end

% IterAIF=AIF_Parker8t(OutAIFParamBeforeFMS,HSampleTs);
% plot(IterAIF,'g','LineWidth',2);
% [[UB; 1]';[LB; 0]';OutAIFParamBeforeFMS;OutAIFParam]
[UB';LB';OutAIFParamBeforeFMS;OutAIFParam]

%%
% gprint(1000,[AIFFinderFN(1:end-4) '.jpg']);
saveas(Handle,[AIFFinderFN(1:end-4) '.fig']);
saveas(Handle,[AIFFinderFN(1:end-4) '.png']);
% close(1000);
Tmp=[AIFFinderFN(1:end-4) '.png'];
AddToLog(WorkingP,'yd_atc','AIF estimation.',Tmp((find(Tmp==filesep,1,'last')+1):end));
%%
HAIF=AIF_Parker9t(OutAIFParam,HSampleTs);
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
% nKeps=100;
% Keps=gpowspace(0,15,nKeps,5)';
% nKeps1=100;
% Keps1I=floor(linspace(1,nKeps,nKeps1));
% 
% HConvd=DCECostFuncgrT1ForConv(HAIF',Keps,HSampleTs,HConvIdxMTriB,HTriB);
% %
% SHConvd=zeros([nTDif nKeps numel(SampleTs)]);
% HSHConvd=zeros([nTDif nKeps numel(HSampleTs)]);
% for i=1:nTDif
%     SHConvd(i,:,:)=interp1(HSampleTs,HConvd',SampleTs+TDif(i),[],'extrap')';
%     HSHConvd(i,:,:)=interp1(HSampleTs,HConvd',HSampleTs+TDif(i),[],'extrap')';
% end
% %
% tic
% [MIdxs CXs Sims] = FindKepBATgAIF(DataToFit,SAIF,SHConvd,Keps1I);
% tFindPK=toc;
% disp(['Finding Kep, BAT took ' num2str(tFindPK)]);

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
%%
WhichToShowAll=randsample(1:nToFit,4*5,true);
% WhichToShowAll=1:20;
for kk=1:5
    figure(Handle+kk);clf;
    set(gcf,'Position',figposition([0 0 100 100]));
    WhichToShow=WhichToShowAll((kk-1)*4+1:kk*4);
    
    for i=1:numel(WhichToShow)
        gsubplot(numel(WhichToShow),i);
        c=WhichToShow(i);
        Regressors=[SAIF(PKs(c,1),:); squeeze(SHConvd2(PKs(c,1),c,:))'];
        HRegressors=[HSAIF(PKs(c,1),:); squeeze(HSHConvd2(PKs(c,1),c,:))'];
        % BAT Kep Vp Ktrans Ve RMS fVal
        CurSim2=((Regressors')*(PKs(c,[3 4])'));
        HCurSim2=((HRegressors')*(PKs(c,[3 4])'));
        ASims(i,:)=CurSim2;
        HASims(i,:)=HCurSim2;
        %         Regressors=[HSAIF(MIdxs(c,1),:); squeeze(HSHConvd(MIdxs(c,1),Keps1I(MIdxs(c,2)),:))'];
        %         CurSim=((Regressors')*(CXs(:,c)));
        %         plot(HSampleTs,CurSim,'b.',SampleTs,CurSim(1:HInterpolationFactor:end),'b*');
        %         plot(SampleTs,Sims(c,:),'r');
        hold on;
        plot(HSampleTs,HCurSim2,'g.');
        plot(SampleTs,CurSim2,'m*');
        plot(SampleTs,DataToFit(c,:),'k-*');
        RError=sqrt(mean((CurSim2'-DataToFit(c,:)).^2))/max(DataToFit(c,:));
        MError=sqrt(PKs(c,7)/nToFit)/max(DataToFit(c,:));
        %             title(['  ' num2str([TShift(c) PKs(c,2) PKs(c,3) PKs(c,4)])]);
        %     title(['  ' num2str([c CostsA(c) Costs(c) WPerCTC(c)])]);
        %         title(PKs(c,7));
        %         title([RError MError]);
        title([num2str(c) ' ' num2str(max(HCurSim2)./max(DataToFit(c,:)),'%2.2f') ' ' num2str(TShift(c))]) ;
        ylabel(max(DataToFit(c,:)));
        set(gca,'YTick',[]);
        set(gca,'XTick',[]);
        %         HSims(i,:)=CurSim2;
        if(i==1)
            xlabel('k:data g:HFim m:Fit');
        end
    end
% end
    %
    %saveas(100+kk,[AIFFinderFN(1:end-4) '_FitOnArt_' num2str(kk) '.fig']);
    %saveas(100+kk,[AIFFinderFN(1:end-4) '_FitOnArt_' num2str(kk) '.png']);
    saveas(Handle+kk,[regexprep(AIFFinderFN(1:end-4),'\.','_') '_FitOnArt_' num2str(kk) '.fig']);
    saveas(Handle+kk,[regexprep(AIFFinderFN(1:end-4),'\.','_') '_FitOnArt_' num2str(kk) '.png']);
%     close(Handle+kk);
%    Tmp=[AIFFinderFN(1:end-4) '_FitOnArt_' num2str(kk) '.png'];
    Tmp=[regexprep(AIFFinderFN(1:end-4),'\.','_') '_FitOnArt_' num2str(kk) '.png'];
    AddToLog(WorkingP,['yd_atd' num2str(kk)],['Fit On chosen voxels example ' num2str(kk)],Tmp((find(Tmp==filesep,1,'last')+1):end));
end
%%
if(nargin>8)
    if(Options.nIterations==0)
        save(AIFFinderFN,'OutAIFParam','ASims','PKs','HASims');
    else
        save(AIFFinderFN,'AllAIF2s','AllAIF2cs','AllRMSs','ARMS','AllAIFPs','OutAIFParam','ASims','PKs');
    end
end
% PKOut=[Keps(Keps1I(MIdxs(:,2))) reshape(TDif(MIdxs(:,1)),[nToFit 1]) CXs(2,:)' CXs(1,:)'];
PKOut=[PKs(:,2) reshape(TDif(PKs(:,1)),[nToFit 1]) PKs(:,4) PKs(:,3)];