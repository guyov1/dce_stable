% Add title to log file
AddToLog(WorkingP,'e_00','\\subsection*{Whole ROI computation.}');

% Struct which holds all the parameters we calculate
TTls={'BATfinal' 'VpFinal' 'KtransFinal' 'Kepfinal' 'VeFinal' 'RSSFinal' 'RSS0' 'RSS1' 'RSS2' 'RSS3' 'F1v0' 'F2v1' 'F3v2' 'BAT1' 'Vp1' 'BAT2' 'Vp2' 'Ktrans2' 'BAT3' 'Vp3' 'Ktrans3' 'Kep3' 'Ve3' 'WhichModel' 'TProb1' 'TProb2' 'TProb3' 'TProbFinal'};

% If using manual chosen artery (without BAT option)
if(Options.MakeNoBATManualArtAnalysis)
    % Add title to log file
    AddToLog(WorkingP,'e_ma','Manual no BAT');
    % Find manual artery
    FindAIFByManualArtNoBAT(WorkingP,true);
end

%% Manual AIF

% If using manual chosen artery (with BAT option and F Min Search for parameters calculation)
if(Options.MakeBATManualArtAnalysis)
    AddToLog(WorkingP,'e_mb','Manual with BAT');
    
    PKMFNman=[WorkingP 'ManualArtBAT_AIFx.mat'];
    PKM3DFNman=[WorkingP 'PKM3DManualx.mat'];
    PKOutP=[WorkingP 'ManualArtBATx' filesep];
    
    load(PKMFNman,'OutAIFParamManual');
    OutAIFParam=OutAIFParamManual;
    
    % Full 3D calculation
    TimeBetweenDCEVolsMin=TimeBetweenDCEVolsFinal/60;
    InterpolationFactor=ceil(TimeBetweenDCEVolsFinal);
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
    %
    clear R1F DCE2D2 adCTC2D B1 Baseline CT1 CosFA DCEM0
    
    % load([WorkingP 'CTCBrainMsk.mat'],'CTC2DBrainMsk');
    % load([WorkingP 'CTCMsk' USStr '.mat'],'CTC2DBigGood','MskCTCGood');
    CTCFN=[WorkingP 'AfterCTC' USStr '.mat'];
    load(CTCFN,'DBrainMask');
    
    MskCTCGood3D=DBrainMask;
    MskCTCGood3D(MskCTCGood3D)=MskCTCGood;
    
    N=size(CTC2DBigGood,1);
    PKs=NaN(N,28);
    NAtATime=5000;
    disp(['There are ' num2str(N) ' voxels to compute ' PKMFNman]);
    before=now;
    
    %     for i=1:NAtATime:N
    %         tic
    %         CurIs=i:min(N,i+NAtATime-1);
    %         PKs(CurIs,:) = FindPKBATgAIFMuraseF4Models_TProb(CTC2DBigGood(CurIs,:),SAIF,SampleTs,CSAIF);
    %         t=toc;
    %         TimeFromStart=now-before;
    %         WillFinishAt=before+TimeFromStart*N/CurIs(end);
    %         disp(['Calculating manual BAT ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
    %     end
    
    % Using parfor for calculation
    num_iter = floor(N/NAtATime);
    PKs_check = zeros(num_iter,NAtATime,size(PKs,2));
    CTC2DBigGood_Sliced = zeros(num_iter,NAtATime,size(CTC2DBigGood,2));
    % Copy CTC2DBigGood into a sliced variable
    for k=1:num_iter
        index = 1 + ( (k-1)*NAtATime );
        CurIs=index:index+NAtATime-1;
        CTC2DBigGood_Sliced(k,:,:) = CTC2DBigGood(CurIs,:);
    end
    
    parfor j=1:num_iter
        tic
        index = 1 + ( (j-1)*NAtATime );
        CurIs=index:index+NAtATime-1;
        
        PKs_check(j,:,:) = FindPKBATgAIFMuraseF4Models_TProb(squeeze(CTC2DBigGood_Sliced(j,:,:)),SAIF,SampleTs,CSAIF);
        %PKs(CurIs,:) = FindPKBATgAIFMuraseF4Models_TProb(CTC2DBigGood(CurIs,:),SAIF,SampleTs,CSAIF);
        
        t=toc;
        TimeFromStart=now-before;
        WillFinishAt=before+TimeFromStart*N/CurIs(end);
        disp(['Calculating manual BAT ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
    end
    
    % Put the resulting paramters in the original variable
    for k=1:num_iter
        index = 1 + ( (k-1)*NAtATime );
        CurIs=index:index+NAtATime-1;
        PKs(CurIs,:) = PKs_check(k,:,:);
    end
    
    % Handle last iteration
    last_index = N - mod(N,NAtATime) + 1;
    tic;
    CurIs=last_index:N;
    PKs(CurIs,:) = FindPKBATgAIFMuraseF4Models_TProb(CTC2DBigGood(CurIs,:),SAIF,SampleTs,CSAIF);
    t=toc;
    TimeFromStart=now-before;
    WillFinishAt=before+TimeFromStart*N/CurIs(end);
    disp(['Calculating manual BAT ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
    
    save(PKM3DFNman,'OutAIFParam','PKs','Options');
    
    % Output to Nii
    mkdir(PKOutP);
    MeanFN=[WorkingP 'DCEMean.nii'];
    disp(PKOutP);
    
    MskX=DBrainMask;
    MskX(MskX)=MskCTCGood;
    
    N=sumn(MskX);
    for i=1:numel(TTls)
        Tmp3D=MskX*0;
        Tmp3D(MskX)=PKs(1:N,i);
        Raw2Nii(Tmp3D,[PKOutP TTls{i} '.nii'],'float32', MeanFN);
    end
    disp('Finished Manual BAT Niftis');
    
    DataP = PKOutP;
    VPstr='Manbat';
    VpNormalization;
    ModelSelections;
end

% If not choosing arteries automatically (default), return.
% Else, execute the code that follows
if(~Options.MakeBATAutoArtAnalysis)
    return;
end

%% Auto AIF

% Get rep voxels as before
% [CVI BinCVI Bin2CVI]=ChooseRepVoxelsForAIFFind(Idx3D>0,CTC2D,BolusStart,Options.Rep_MaxAroundBolus,Options.Rep_RatioToEnd,Options.Rep_nPerSet);
% [CVI BinCVI Bin2CVI]=ChooseRepVoxelsForAIFFind(Idx3D>0,CTC2D,BolusStart,5,5,2);
%
% % Options=struct('SubSecRes',{[]},'MaxTDif',{[]});

% Get the noise of each voxel
% DataNoise=rmadCTC2D(CVI);
% Get the date to fit to
% DataToFit=CTC2D(CVI,:);

% Get the optimal paramters computed before
PKMFN   = [WorkingP 'PKM' USStr '.mat'];
PKM3DFN = [WorkingP 'PKM3D' USStr '.mat'];
load(PKMFN,'OutAIFParam');

InspectedAIFParamsFN=[WorkingP 'InspectedAIFParams.mat'];
if(exist(InspectedAIFParamsFN,'file'))
    AddToLog(WorkingP,'e_00ins','Using inspected AIF Parameters.');
    Tmp=load(InspectedAIFParamsFN);
    OutAIFParam=Tmp.InspectedParams;
end

%-------------- Full 3D calculation ----------------

GoodTs=Sts;
GoodTIdxs=1:numel(Sts);
if (Num_Of_T1_Maps>1)
    NumBefore=size(Additional_before_main,2);
    
    % If SingleAngleIdxs was not defined in DCET1)CTCf, set it as '0'
    if ( ~exist('SingleAngleIdxs') )
        SingleAngleIdxs = 0;
    end
    
    GoodTIdxs=[(NumBefore+1):(numel(ts)+NumBefore) numel(ts)+SingleAngleIdxs-1];
    GoodTs=[Sts(1:numel(ts)) Additional_T1_Maps_Time_Diff_Sets(SingleAngleIdxs)'*86400/60];
    
    GoodTIdxs=1:numel(SampleTsNoBefore);
    GoodTs=SampleTsNoBefore/60;
    % figure;plot(GoodTs,DataToFit(1:10:end,GoodTIdxs)','.-');
end
SampleTs=GoodTs;
TimeVec=SampleTs;
% As before, set all relevant paramters
TimeBetweenDCEVolsMin=TimeBetweenDCEVolsFinal/60;
InterpolationFactor=ceil(TimeBetweenDCEVolsFinal);

HInterpolationFactor = ceil(InterpolationFactor*Options.SubSecondResolution);
% HInterpolationFactor=ceil(InterpolationFactor*Options.SubSecRes);
Hdt                  = TimeBetweenDCEVolsMin/HInterpolationFactor;
HSampleTs            = 0:Hdt:SampleTs(end);
DTimeVec=diff(TimeVec);
F=find(DTimeVec>DTimeVec(1)*2,1);
nNormalTs=numel(TimeVec);
if(~isempty(F))
    nNormalTs=F(1);
end
HSampleTs=[0:Hdt:SampleTs(nNormalTs) TimeVec((nNormalTs+1):end)];

ThreeSec = ceil(Options.MaxTDif_ForWholeVOI/(Hdt*60));
TDif     = -Hdt*ThreeSec:Hdt:Hdt*ThreeSec;
nTDif    = numel(TDif);

AddToLog(WorkingP,'e_t1',['HInterpolationFactor: ' num2str(HInterpolationFactor)]);

% AIF_Parker8t=@(x,t) AIF_Parkerg2( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8))*x(2);

InspectedAIFParamsTimeFN=[WorkingP 'InspectedAIFParamsTime.mat'];
if(exist(InspectedAIFParamsTimeFN,'file'))
    
    Tmp1=load(InspectedAIFParamsTimeFN);
    OldTimeVec=Tmp1.InspectedParamsTimeSamples;
    AIF_Parker9t=@(x,t) AIF_Parkerg3( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8),x(9),max(OldTimeVec))*x(2);
else
    AIF_Parker9t=@(x,t) AIF_Parkerg3( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8),x(9),max(TimeVec))*x(2);
end



CSVOut=[AIF_Parker9t(OutAIFParam,SampleTs); SampleTs];
CSVOutFN=[WorkingP 'AIF.csv'];
csvwrite(CSVOutFN,CSVOut);
HTForCSV=0:0.1:max(ceil(SampleTs));
HCSVOut=[AIF_Parker9t(OutAIFParam,HTForCSV); HTForCSV];
HCSVOutFN=[WorkingP 'HTR-AIF.csv'];
csvwrite(HCSVOutFN,HCSVOut);
% Create Parker's AIF
HAIF=AIF_Parker9t(OutAIFParam,HSampleTs);
% Create Cummulative Parker's AIF
CHAIF = cumtrapz(HSampleTs,HAIF);

SAIF  = zeros([nTDif numel(SampleTs)]);
CSAIF = zeros([nTDif numel(SampleTs)]);

for i=1:nTDif
    SAIF(i,:)  = interp1(HSampleTs,HAIF,SampleTs+TDif(i),[],'extrap');
    CSAIF(i,:) = interp1(HSampleTs,CHAIF,SampleTs+TDif(i),[],'extrap');
end

% Clear unnecessary old variables
clear R1F DCE2D2 adCTC2D B1 Baseline CT1 CosFA DCEM0

% load([WorkingP 'CTCBrainMsk.mat'],'CTC2DBrainMsk');
%load([WorkingP 'CTCMsk' USStr '.mat'],'CTC2DBigGood','MskCTCGood');

CTCFN = [WorkingP 'AfterCTC' USStr '.mat'];
load(CTCFN,'DBrainMask');

CTC2DBigGood=CTC2DBigGood(:,GoodTIdxs);

MskCTCGood3D=DBrainMask;
MskCTCGood3D(MskCTCGood3D)=MskCTCGood;

N        = size(CTC2DBigGood,1);
PKs      = NaN(N,28);
NAtATime = 5000; % Size of each chunck of voxels to work on
disp(['There are ' num2str(N) ' voxels to compute ' PKMFN]);
before   = now;

% Go over all voxels and create for each parameters using AIF
% Use parfor calculation for better performance
num_iter            = floor(N/NAtATime);
PKs_check           = zeros(num_iter,NAtATime,size(PKs,2));
CTC2DBigGood_Sliced = zeros(num_iter,NAtATime,size(CTC2DBigGood,2));

% Copy CTC2DBigGood into a sliced variable
for k=1:num_iter
    index = 1 + ( (k-1)*NAtATime );
    CurIs=index:index+NAtATime-1;
    CTC2DBigGood_Sliced(k,:,:) = CTC2DBigGood(CurIs,:);
end
%%
% Executing parameters calculation
parfor j=1:num_iter
    tic;
    index = 1 + ( (j-1)*NAtATime );
    CurIs = index:index+NAtATime-1;
    
    % Calculate parameters for the current chunck of voxels
    PKs_check(j,:,:) = FindPKBATgAIFMuraseF4Models_TProb(squeeze(CTC2DBigGood_Sliced(j,:,:)),SAIF,SampleTs,CSAIF);
    
    t=toc;
    TimeFromStart=now-before;
    WillFinishAt=before+TimeFromStart*N/CurIs(end);
    disp(['Calculating auto AIF BAT ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
end
%% Put the resulting paramters in the original variable
for k=1:num_iter
    index        = 1 + ( (k-1)*NAtATime );
    CurIs        = index:index+NAtATime-1;
    PKs(CurIs,:) = PKs_check(k,:,:);
end

% Handle last iteration
last_index   = N - mod(N,NAtATime) + 1;
tic;
CurIs        = last_index:N;
% Calculate parameters for the current chunck of voxels
PKs(CurIs,:) = FindPKBATgAIFMuraseF4Models_TProb(CTC2DBigGood(CurIs,:),SAIF,SampleTs,CSAIF);
t = toc;
TimeFromStart = now - before;
WillFinishAt  = before+TimeFromStart*N/CurIs(end);
disp(['Calculating auto AIF BAT ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
% Save calculated parameters
save(PKM3DFN,'OutAIFParam','PKs');

% Output to Nii
PKOutP = [WorkingP 'AutoArtBAT' filesep];
mkdir(PKOutP);
MeanFN = [WorkingP 'DCEMean.nii'];
disp(PKOutP);

MskX       = DBrainMask;
MskX(MskX) = MskCTCGood;

N=sumn(MskX);
for i=1:numel(TTls)
    Tmp3D       = MskX*0;
    Tmp3D(MskX) = PKs(1:N,i);
    Raw2Nii(Tmp3D,[PKOutP TTls{i} '.nii'],'float32', MeanFN);
end
disp('Finished Niftis');

%% Normalize Vp values (ASK GILAD HOW IT WORKS)
DataP = PKOutP;
VPstr='Autobat';
VpNormalization;
ModelSelections;
BatFinalFN=[DataP 'BATfinal.nii'];
Tmp=loadniidata(BatFinalFN);
TDifx=[-100 -TDif];
Tmp(isfinite(Tmp))=TDifx(Tmp(isfinite(Tmp))+1)*60;
BatFinalFNsec=[DataP 'BATfinalSec.nii'];
Raw2Nii(Tmp,BatFinalFNsec,'float32', MeanFN);
[XX, PeakInd]=max(HAIF);
AIFPeakTime=HSampleTs(PeakInd)*60;
BatFinalFNsecN=[DataP 'BATfinalSecN.nii'];
Raw2Nii(AIFPeakTime+Tmp,BatFinalFNsecN,'float32', MeanFN);