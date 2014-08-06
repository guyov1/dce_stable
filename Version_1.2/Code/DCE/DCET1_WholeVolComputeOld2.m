if(Options.MakeNoBATManualArtAnalysis)
    FindAIFByManualArtNoBAT(WorkingP)
end

if(Options.MakeBATManualArtAnalysis)
    PKMFNman=[WorkingP 'ManualArtBAT_AIF.mat'];
    PKM3DFNman=[WorkingP 'PKM3DManual.mat'];

    load(PKMFNman,'OutAIFParamManual');
    OutAIFParam=OutAIFParamManual;
    %% Full 3D calculation
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
    clear R1F DCE2D2 adCTC2D B1 Baseline CT1 CosFA DCEM0

    % load([WorkingP 'CTCBrainMsk.mat'],'CTC2DBrainMsk');
    load([WorkingP 'CTCMsk' USStr '.mat'],'CTC2DBigGood','MskCTCGood');
    CTCFN=[WorkingP 'AfterCTC' USStr '.mat'];
    load(CTCFN,'DBrainMask');

    MskCTCGood3D=DBrainMask;
    MskCTCGood3D(MskCTCGood3D)=MskCTCGood;

    N=size(CTC2DBigGood,1);
    PKs=NaN(N,9);
    NAtATime=5000;
    disp(['There are ' num2str(N) ' voxels to compute ' PKMFNman]);
    before=now;
    for i=1:NAtATime:N
        tic
        CurIs=i:min(N,i+NAtATime-1);
        PKs(CurIs,:) = FindPKBATgAIFMuraseF_TProb(CTC2DBigGood(CurIs,:),SAIF,SampleTs,CSAIF);
        t=toc;
        TimeFromStart=now-before;
        WillFinishAt=before+TimeFromStart*N/CurIs(end);
        disp(['Calculating manual ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
    end
    save(PKM3DFNman,'OutAIFParam','PKs');
    %% Output to Nii
    MeanVol=loadniidata(MeanFN);

    % BAT Kep Vp Ktrans Ve RMS fVal
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
    RMS3D(MskCTCGood3D)=PKs(:,6);
    F3D=MskCTCGood3D*0;
    F3D(MskCTCGood3D)=PKs(:,7);
    F03D=MskCTCGood3D*0;
    F03D(MskCTCGood3D)=PKs(:,8);

    TrgP=[WorkingP 'ManualArtBAT' filesep];
    mkdir(TrgP);
    Raw2Nii(Kep3D,[TrgP 'Kep' USStr '.nii'],'float32', MeanFN);
    Raw2Nii(BAT3D,[TrgP 'BAT' USStr '.nii'],'float32', MeanFN);
    Raw2Nii(Vp3D,[TrgP 'Vp' USStr '.nii'],'float32', MeanFN);
    Raw2Nii(Ktrans3D,[TrgP 'Ktrans' USStr '.nii'],'float32', MeanFN);
    Raw2Nii(Ve3D,[TrgP 'Ve' USStr '.nii'],'float32', MeanFN);
    Raw2Nii(RMS3D,[TrgP 'RMS' USStr '.nii'],'float32', MeanFN);
    Raw2Nii(F3D,[TrgP 'F' USStr '.nii'],'float32', MeanFN);
    Raw2Nii(F03D,[TrgP 'F0' USStr '.nii'],'float32', MeanFN);
end

if(~Options.MakeBATAutoArtAnalysis)
    return;
end

% [CVI BinCVI Bin2CVI]=ChooseRepVoxelsForAIFFind(Idx3D>0,CTC2D,BolusStart,Options.Rep_MaxAroundBolus,Options.Rep_RatioToEnd,Options.Rep_nPerSet);
% % [CVI BinCVI Bin2CVI]=ChooseRepVoxelsForAIFFind(Idx3D>0,CTC2D,BolusStart,5,5,2);
%
% % Options=struct('SubSecRes',{[]},'MaxTDif',{[]});
%
% DataNoise=rmadCTC2D(CVI);
% DataToFit=CTC2D(CVI,:);

PKMFN=[WorkingP 'PKM' USStr '.mat'];
PKM3DFN=[WorkingP 'PKM3D' USStr '.mat'];

%%
load(PKMFN,'OutAIFParam');
%% Full 3D calculation
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
clear R1F DCE2D2 adCTC2D B1 Baseline CT1 CosFA DCEM0

% load([WorkingP 'CTCBrainMsk.mat'],'CTC2DBrainMsk');
load([WorkingP 'CTCMsk' USStr '.mat'],'CTC2DBigGood','MskCTCGood');
CTCFN=[WorkingP 'AfterCTC' USStr '.mat'];
load(CTCFN,'DBrainMask');

MskCTCGood3D=DBrainMask;
MskCTCGood3D(MskCTCGood3D)=MskCTCGood;

N=size(CTC2DBigGood,1);
PKs=NaN(N,9);
NAtATime=5000;
disp(['There are ' num2str(N) ' voxels to compute ' PKMFN]);
before=now;
for i=1:NAtATime:N
    tic
    CurIs=i:min(N,i+NAtATime-1);
    PKs(CurIs,:) = FindPKBATgAIFMuraseF_TProb(CTC2DBigGood(CurIs,:),SAIF,SampleTs,CSAIF);
    t=toc;
    TimeFromStart=now-before;
    WillFinishAt=before+TimeFromStart*N/CurIs(end);
    disp(['Calculating ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
end
save(PKM3DFN,'OutAIFParam','PKs');
%% Output to Nii
MeanVol=loadniidata(MeanFN);

% BAT Kep Vp Ktrans Ve RMS fVal
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
RMS3D(MskCTCGood3D)=PKs(:,6);
F3D=MskCTCGood3D*0;
F3D(MskCTCGood3D)=PKs(:,7);
F03D=MskCTCGood3D*0;
F03D(MskCTCGood3D)=PKs(:,8);

Raw2Nii(Kep3D,[WorkingP 'Kep' USStr '.nii'],'float32', MeanFN);
Raw2Nii(BAT3D,[WorkingP 'BAT' USStr '.nii'],'float32', MeanFN);
Raw2Nii(Vp3D,[WorkingP 'Vp' USStr '.nii'],'float32', MeanFN);
Raw2Nii(Ktrans3D,[WorkingP 'Ktrans' USStr '.nii'],'float32', MeanFN);
Raw2Nii(Ve3D,[WorkingP 'Ve' USStr '.nii'],'float32', MeanFN);
Raw2Nii(RMS3D,[WorkingP 'RMS' USStr '.nii'],'float32', MeanFN);
Raw2Nii(F3D,[WorkingP 'F' USStr '.nii'],'float32', MeanFN);
Raw2Nii(F03D,[WorkingP 'F0' USStr '.nii'],'float32', MeanFN);