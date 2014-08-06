if(Options.MakeNoBATManualArtAnalysis)
    FindAIFByManualArtNoBAT(WorkingP)
end

if(Options.MakeBATManualArtAnalysis)
    PKMFNman=[WorkingP 'ManualArtBAT_AIF.mat'];
    PKM3DFN=[WorkingP 'PKM3DManual.mat'];

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
    PKs=NaN(N,24);
    NAtATime=5000;
    disp(['There are ' num2str(N) ' voxels to compute ' PKMFNman]);
    before=now;
    for i=1:NAtATime:N
        tic
        CurIs=i:min(N,i+NAtATime-1);
%         PKs(CurIs,:) = FindPKBATgAIFMuraseF_TProb(CTC2DBigGood(CurIs,:),SAIF,SampleTs,CSAIF);
        PKs(CurIs,:) = FindPKBATgAIFMuraseF4Models_TProb(CTC2DBigGood(CurIs,:),SAIF,SampleTs,CSAIF);
        t=toc;
        TimeFromStart=now-before;
        WillFinishAt=before+TimeFromStart*N/CurIs(end);
        disp(['Calculating manual ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
    end
    save(PKM3DFN,'OutAIFParam','PKs');
    disp(['Finished manual ' WorkingP]);
    
    % Output to Nii
    TTls={'BATfinal' 'VpFinal' 'KtransFinal' 'Kepfinal' 'VeFinal' 'RSSFinal' 'RSS0' 'RSS1' 'RSS2' 'RSS3' 'F1v0' 'F2v1' 'F3v2' 'BAT1' 'Vp1' 'BAT2' 'Vp2' 'Ktrans2' 'BAT3' 'Vp3' 'Ktrans3' 'Kep3' 'Ve3' 'WhichModel'};

    RelaxP=[WorkingP 'ManualArtBAT' filesep];
    mkdir(RelaxP);
    MeanFN=[WorkingP 'DCEMean.nii'];
    disp(RelaxP);
    
    MskX=DBrainMask;
    MskX(MskX)=MskCTCGood;

    N=sumn(MskX);
    for i=1:numel(TTls)
        Tmp3D=MskX*0;
        Tmp3D(MskX)=PKs(1:N,i);
        Raw2Nii(Tmp3D,[RelaxP TTls{i} '.nii'],'float32', MeanFN);
    end
    disp(['Finished nii ' WorkingP]);
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
PKs=NaN(N,24);
NAtATime=5000;
disp(['There are ' num2str(N) ' voxels to compute ' PKMFN]);
before=now;
for i=1:NAtATime:N
    tic
    CurIs=i:min(N,i+NAtATime-1);
    PKs(CurIs,:) = FindPKBATgAIFMuraseF4Models_TProb(CTC2DBigGood(CurIs,:),SAIF,SampleTs,CSAIF);
    t=toc;
    TimeFromStart=now-before;
    WillFinishAt=before+TimeFromStart*N/CurIs(end);
    disp(['Calculating ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
end
save(PKM3DFN,'OutAIFParam','PKs');

%% Output to Nii
TTls={'BATfinal' 'VpFinal' 'KtransFinal' 'Kepfinal' 'VeFinal' 'RSSFinal' 'RSS0' 'RSS1' 'RSS2' 'RSS3' 'F1v0' 'F2v1' 'F3v2' 'BAT1' 'Vp1' 'BAT2' 'Vp2' 'Ktrans2' 'BAT3' 'Vp3' 'Ktrans3' 'Kep3' 'Ve3' 'WhichModel'};

RelaxP=[WorkingP 'ManualArtBAT' filesep];
mkdir(RelaxP);
MeanFN=[WorkingP 'DCEMean.nii'];
disp(RelaxP);

MskX=DBrainMask;
MskX(MskX)=MskCTCGood;

N=sumn(MskX);
for i=1:numel(TTls)
    Tmp3D=MskX*0;
    Tmp3D(MskX)=PKs(1:N,i);
    Raw2Nii(Tmp3D,[RelaxP TTls{i} '.nii'],'float32', MeanFN);
end
disp(['Finished nii ' WorkingP]);
