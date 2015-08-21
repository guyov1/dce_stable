
function DCET1_Prepare4Df(CurMainDCEInfo,WorkingP,Force, Options)

% Fix Path according to windows/unix 
if ( strcmp(filesep,'\') ) % Windows
    CurMainDCEInfo.Path=strrep(CurMainDCEInfo.Path,'\\','\');    
else % Unix
    CurMainDCEInfo.Path=strrep(CurMainDCEInfo.Path,'//','/');        
end

% Get DCE Main Path
DCEMainP=CurMainDCEInfo.Path;

% Get the main DICOM data
CurFullInfo=dicominfo(CurMainDCEInfo.Filename);
Philips=strhas(CurMainDCEInfo.Filename,'Uptake');
%% Check if stage was done already / "force" was used.

% Path for the output data calculated in this stage
PrepareFN=[WorkingP 'AfterPrepare4D.mat'];

% If the data already exists, don't calculate again
% if(exist(PrepareFN,'file') && ~Force)
%     D=dir(PrepareFN);
% %     if(D.datenum>DDate)
% %         loadBut(PrepareFN,{CurVars.name});
%         disp('DCE_Prepeare4D already computed');
%         return;
% %     end
% end

% If the file exists and we got here, we used Force -> Delete anyway.
if(exist(PrepareFN,'file'))
%     delete(PrepareFN);
end
LogFN=[WorkingP 'Log.mat'];
delete(LogFN)
SN=[ToShortName(CurMainDCEInfo.Name) ' ' CurMainDCEInfo.SeriesDate];
Log.a_00={['\\title{' SN '}\r\n\\maketitle\r\n']};
save(LogFN,'Log');

AddToLog(WorkingP,'a_01','\\subsection*{Preprocess}');
AddToLog(WorkingP,'z_aaa', '\\newpage\r\n\\subsubsection*{Options}');
for i=fieldnames(Options)'
    AddToLog(WorkingP,['zz_a' i{1}],strrep([i{1} ' = ' num2str(Options.(i{1}))],'_',':'));
end
AddToLog(WorkingP,'zz_b','---------------------');

%TimeBetweenDCEVolsPerSlice=double(min(CurFullInfo.AcquisitionMatrix(CurFullInfo.AcquisitionMatrix>0))*CurFullInfo.RepetitionTime)/(1000*max(CurFullInfo.EchoTrainLength,1));
AddToLog(WorkingP,'zza_2a0',['AcquisitionMatrix: ' num2str(CurFullInfo.AcquisitionMatrix(CurFullInfo.AcquisitionMatrix>0)')]);
AddToLog(WorkingP,'zza_2a1',['EchoTrainLength: ' num2str(max(CurFullInfo.EchoTrainLength,1))]);
AddToLog(WorkingP,'zza_2a2',['RepetitionTime: ' num2str(CurFullInfo.RepetitionTime)]);
AddToLog(WorkingP,'zza_2a3',['PercentPhaseFieldOfView: ' num2str(CurFullInfo.PercentPhaseFieldOfView)]);
AddToLog(WorkingP,'zza_2a4',['EchoTime: ' num2str(CurFullInfo.EchoTime)]);
AddToLog(WorkingP,'zza_2a5',['FlipAngle: ' num2str(CurFullInfo.FlipAngle)]);

% Remove all previous directories
if(exist([WorkingP 'Seg_DCEMean'],'file'))
    rmdir([WorkingP 'Seg_DCEMean'],'s');
end
if(exist([WorkingP 'Seg_DCEMeanL'],'file'))
    rmdir([WorkingP 'Seg_DCEMeanL'],'s');
end
if(exist([WorkingP 'DCEMainCoregedLongitudinal'],'file'))
    rmdir([WorkingP 'DCEMainCoregedLongitudinal'],'s');
end

% Display the starting time and date of the process
disp(['DCE_Prepeare4D, start ' datestr(now)]);

%% Initial parameters values

% ASK GILAD - what is the purpose of the following parameters
% Asnwer - Min singal from each the data is relevant.
%          Min enhancment for mask usage (less than that is not interesting)
MinSignal=50;
MinEnhancementR=1.02;
% MinSPMBrainValue=0.3;
%% Step 1 - Clean and realign DCE images

disp('Dicom to Nifti');
DCEMNiiOrigP=[WorkingP 'DCEMainNii' filesep];
mkdir(DCEMNiiOrigP);

% Get all the nifty files
DDCE=dir([DCEMNiiOrigP '*.nii']);

if(isempty(DDCE) || Force)
    % Convert Dicom to Nifty
    gDicom2Nifti(DCEMainP,[DCEMNiiOrigP 'vol.nii']);
    DDCE=dir([DCEMNiiOrigP '*.nii']);
    if(isempty(DDCE))
        error('DCET1_Prepare4Df - Dicom2Nifi didnt work');
    end
end

% ASK GILAD
% Clean - not needed?

% Realign and extract mean
disp('Coregistration');
% The following coregister the images.
% ASK GILAD - did he write this function? if so, explain the steps.
DCE_Coreg;
DCEMNiiP=DCECoregP;
% if(CoregLongPerformed)
%     DCE_CoregLong;
% end
%% Step 2 - Load all and create 4D and mean files

disp('Load 4D and calculate mean');
DDCE=dir([DCEMNiiP '*.nii']);
DCEFNs=strcat(DCEMNiiP,{DDCE.name})';
DCE4D=loadCoregedNiftis(DCEFNs(1:end-Philips));
SDCE=size(DCE4D);
nVols=size(DCE4D,4);
DCE4DFN=[WorkingP 'DCE4D.nii'];
Raw2Nii(DCE4D,DCE4DFN,'uint16',DCEFNs{1});
MeanVol=mean(DCE4D,4);
MeanFN=[WorkingP 'DCEMean.nii'];

% if(CoregLongPerformed)
%     MeanFN=[WorkingP 'DCEMeanL.nii'];
% end
% copyfile([DCEMNiiP 'vol_' num2str(floor(size(DCE4D,4)/2)) '.nii'],MeanFN,'f');
Raw2Nii(MeanVol,MeanFN,'float32',DCEFNs{1});

% ASK GILAD - what is the purpose of the following?
% ASNWER - the time between each slice is dependent on the resolution, repetition time and acceleration ingredients.
if(Philips)
    TimeBetweenDCEVols=15;
else
    if(isfield(CurFullInfo,'NumberOfTemporalPositions') && isfield(CurFullInfo,'Private_0019_105a'))
        if(numel(CurFullInfo.Private_0019_105a)==1)
            A=CurFullInfo.Private_0019_105a;
        else
            Int8x4ToFloat;
        end
        TimeBetweenDCEVols=(A./[CurFullInfo.NumberOfTemporalPositions])/1e6;
        AddToLog(WorkingP,'a_2b3','Used new time calculation');
    else
        TimeBetweenDCEVolsPerSlice=double(min(CurFullInfo.AcquisitionMatrix(CurFullInfo.AcquisitionMatrix>0))*CurFullInfo.RepetitionTime)/(1000*max(CurFullInfo.EchoTrainLength,1))*CurFullInfo.PercentPhaseFieldOfView/100.;
        TimeBetweenDCEVols=TimeBetweenDCEVolsPerSlice*SDCE(3);
        AddToLog(WorkingP,'a_2b3','Used old time calculation');
        AddToLog(WorkingP,'a_2b3b',['Old time calc: ' num2str(TimeBetweenDCEVols)]');
    end
end

if(isfield(Options,'TimeMultiplier'))
    TimeBetweenDCEVols=TimeBetweenDCEVols*Options.TimeMultiplier;
end
if(isfield(Options,'TimeMultiplier') && Options.TimeMultiplier<0)
    TimeBetweenDCEVols=-Options.TimeMultiplier;
end
TimeBetweenDCEVolsMin=TimeBetweenDCEVols/60;
%% Step 3 - Create mutual information curve to see when the registration did not succeed

disp('Calculate mutual information');
% Mutual information curve
% ASK GILAD - explain how he calculated the MI and what it means exactly.
% MIC=CalcMICurve(DCEFNs,MeanFN);

%%  Step 4 - Show movements -> 	Create a montage of all middle slices during time.

disp('Display coreged');
% if(false)
    % Middle Slice
    MidSli=floor(SDCE(3)/2);
    
    % Dynamic Range?
    % ASK GILAD - why does he divide by 3?
    DR=[min(DCE4D(:)) max(DCE4D(:))/3];
        figure(987);clf;
        
    % Cretae a montage image of all the middle slices during time.
    montage(mritransformNoSqueeze((DCE4D(:,:,MidSli,:))-DR(1))./(DR(2)-DR(1)));
    saveas(987,[WorkingP 'CoregedMidSlice.png']);
    saveas(987,[WorkingP 'CoregedMidSlice.fig']);
    close(987);
% end
AddToLog(WorkingP,'xa_1','Coreged.','CoregedMidSlice.png');
%% Step 5 -  Masks creation

% TimeBetweenDCEVols
% TimeBetweenDCEVols=TimeBetweenDCEVolsPerSlice*SDCE(3);
% TimeBetweenDCEVolsMin=TimeBetweenDCEVols/60;
AddToLog(WorkingP,'a_2a',['Size:' num2str(SDCE)]);
AddToLog(WorkingP,'a_2b',['Time between DCE vols:' num2str(TimeBetweenDCEVols')]);
if(abs(TimeBetweenDCEVols-6)>1)
    AddToLog(WorkingP,'a_2b1','Warning - weird time between DCE vols.',[],2);
end
% Frequency of images over time
% The following variable is used later in DCET1_CTCf.m
Min2Timestep=1/TimeBetweenDCEVolsMin;

%% Min signal mask

%  Mask all time volumes to see whether in the entire test, the signal was stronger than the minimal.
MskMinSignal=all(DCE4D>MinSignal,4);

%% Masking by FSL (previously used SPM)
disp('Brain extraction');

% Running FSL's bet function to extract the brain out of the image (mask)
if(Options.Mask_thresh<0)
    [BetStrippedMeanFN, BetMaskFN]=mcbet2(MeanFN,Force,abs(Options.Mask_thresh));
else
    [BetStrippedMeanFN, BetMaskFN]=mcbet2(MeanFN,Force);
end
ManMaskFN=[WorkingP 'Manual_BrainMask.nii'];
if(exist(ManMaskFN,'file'))
    %AddToLog(WorkingP,'a_2caaaaa','Using Manual_BrainMask');
    AddToLog(WorkingP,'a_2caaaaa','Using Manual.BrainMask');
    BrainMask=loadniidata(ManMaskFN)>0;
else
    BrainMask=loadniidata(BetMaskFN)>0;
end


%% Masking bad slices

%  Mark the first and the last slices to be false
BadSlicesA=[1 size(MeanVol,3)];

% Remove the first and last slices
BrainMask(:,:,BadSlicesA)=false;

%Selecting slices
disp('Find problematic slices');

MedSlice=zeros(SDCE(3),1);
% Go over every slice
for i=1:SDCE(3)
    % Slice mask at every iteration holds values different than 0 for the
    % current slice only (this is why we initiate it with zeros)
    SliceMsk=zeros(SDCE(1),SDCE(2),SDCE(3));
    SliceMsk(:,:,i)=BrainMask(:,:,i);
    
    % Take the first time index only because we want an image before the enhancement.
    % We should have taken the entire baseline (all images before
    % enhancmenet), but this is an approximation
    Tmp=Reshape4d22d(DCE4D(:,:,:,1),SliceMsk);
    % Get the median value of every slice
    MedSlice(i,:)=median(Tmp,1);
end

% The following function is problematic (dimenstions for repPlus,repMulti
% are incorrect)
%NMedSlice=repMulti(repPlus(MedSlice,-mean(MedSlice,1)),1./std(MedSlice,0,1));

% BadSlices=any(abs(NMedSlice)>2,2);
NaNSlices=isnan(MedSlice);
disp('Looking for weird slices');
if(numel(MedSlice(~NaNSlices))<2)
    error('Only one slice!');
end

BadSlices=NaNSlices;

% The following function is checking if there are any other bad slices
% A slice is considered bad if its median is far (>50) than the 3 middle slices (which should have a reliable signal
BadSlices(~NaNSlices)=abs(MedSlice(~NaNSlices)-mean(MedSlice(MidSli-1:MidSli+1)))>50-Philips*30;
BadSlicesF=find(BadSlices);
disp(['Ignoring slices ' num2str(BadSlicesF')]);
AddToLog(WorkingP,'a_2ccc',['Ignoring slices ' num2str(BadSlicesF')]);
figure(78362);clf;subplot(1,2,1);
% Plot the median value of each slice
plot(1:SDCE(3),MedSlice,'b',BadSlicesF,MedSlice(BadSlicesF),'ro');
title('Median of slices and ignored ones');
xlabel('Slices');
ylabel('Median');
MskSlices=ones(size(MskMinSignal))>0;
MskSlices(:,:,BadSlicesF)=false;
% Get a list of all bad slices (unite the extreme slices and the ones >50)
BadSlicesF2=union(BadSlicesF,[1 size(MeanVol,3)]);
BrainMskFN=[WorkingP 'BrainMask.nii'];
BrainMask(:,:,BadSlicesF2)=0;
Raw2Nii(BrainMask,BrainMskFN,'float32',MeanFN);
MeanVolGoodSlices=MeanVol;
MeanVolGoodSlices(~BrainMask)=0;
MeanVolGoodSlices(:,:,BadSlicesF2)=0;
MeanFNGoodSlices=[WorkingP 'MeanFNNo' GroupToStr(num2strC(BadSlicesF2)) '.nii'];
Raw2Nii(MeanVolGoodSlices,MeanFNGoodSlices,'float32', MeanFN);
%% Step 6 - Find bolus start, compute baseline

disp('Rough estimation of bolus time');
% We mask again for all slices with value bigger than minimum (this time for all time periods and not just the first)
DCE2D=Reshape4d22d(DCE4D,MskMinSignal);
% Get the median of each of the 3d images for the entire time slots 
MedTC=median(DCE2D,1);
% [~, optmixture] = GaussianMixture(MedTC', 3, 0,false);
% [Z Grouped]=max(optmixture.pnk,[],2);
% if(max(Grouped)==1)
%     Smoothed=conv2(MedTC,ones(1,ceil(nVols/10)),'same');
%     [a BolusStart]=max(Smoothed);
% else
%     Grouped(end)=Grouped(1)+1;
%     BolusStart=find(Grouped~=Grouped(1),1);end
% Second method
TwoMinTimePoint=floor(2/TimeBetweenDCEVolsMin);
Ps=zeros(1,numel(MedTC))+2;
% We use the t-test to get the biggest probability that the distribution of the sample
% is diffrent than the rest of the test ( -> smallest Ps value)
for i=3:min(TwoMinTimePoint,numel(MedTC)-2) %Take the minimum out of 2 minutes frame to 2 frames before the end
    [h Ps(i)]=ttest2(MedTC(1:i),MedTC((i+1):end),[],[],'unequal');
end
mLPs=-log(Ps);
% figure;plot(1:numel(MedTC),MedTC,'b',1:numel(MedTC),mLPs.*(max(MedTC)-min(MedTC))./(max(mLPs)-min(mLPs))+min(MedTC),'r')
[Tmp, BolusStart]=max(mLPs);

% ASK GILAD - why did he add +1?
BolusStart=BolusStart+1;
BolusStartMin=(BolusStart-1)*TimeBetweenDCEVolsMin;
% BolusStart=find(MedTC>MedTC(1)+20,1);
% The base line is the mean of the first images until the bolus
Baseline=mean(DCE4D(:,:,:,1:(BolusStart-2)),4);
BaselineFN=[WorkingP 'Baseline.nii'];
Raw2Nii(Baseline,BaselineFN,'float32', MeanFN);
figure(78362);subplot(1,2,2);
plot(MedTC); hold on;plot([BolusStart BolusStart],[min(MedTC) max(MedTC)],'r');
title('Bolus start approximation');

%% Step 6.3 new SPM masking with quick B1
DCEMeanSegQB1=SPM_SegmentWithB1(MeanFN,Force);
FMaskFN=[WorkingP 'FBrainMsk.nii'];

% If the brain mask have black holes in it, fill them
FBrainMask=bwfillHoles3Dby2D(BrainMask);
A=load_untouch_nii(DCEFNs{1});
A.img=int16(FBrainMask);
save_untouch_nii(A,FMaskFN);
disp('FBrainMsk finished');
%% Step 6.5 SPM masking
BaselineNoBadSliFN=[WorkingP 'BaselineNoBadSli.nii'];
BaselineNoBadSli=Baseline;
BaselineNoBadSli(:,:,BadSlicesF2)=NaN;
Raw2Nii(BaselineNoBadSli,BaselineNoBadSliFN,'float32', MeanFN);
% 
% DCEMeanSegP=SPM_Segment(BaselineNoBadSliFN,Force,[],false);
DCEMeanSegP=DCEMeanSegQB1;
C1=loadniidata([DCEMeanSegP 'c1ForSeg.nii'])/256;
C2=loadniidata([DCEMeanSegP 'c2ForSeg.nii'])/256;
C3=loadniidata([DCEMeanSegP 'c3ForSeg.nii'])/256;
if(Options.Mask_thresh>0)
    MinSPMBrainValue=Options.Mask_thresh;
else
    MinSPMBrainValue=0.5;    
end
% % BrainMask=(C1+C2)>MinSPMBrainValue;
% % 
BrainMaskA=(C1+C2+C3)>MinSPMBrainValue;
BrainMaskA=bwfillHoles3Dby2D(BrainMaskA);
se=strel('disk',2,8);
BrainMaskA=imerode(BrainMaskA,se);
BrainMaskA=bwfillHoles3Dby2D(BrainMaskA);
if((Options.Mask_thresh>0) && ~exist(ManMaskFN,'file'))
    BrainMask=BrainMaskA;
    FBrainMask=BrainMask;
end
Raw2Nii(BrainMask,BrainMskFN,'float32',MeanFN);
Raw2Nii(BrainMask,FMaskFN,'float32',MeanFN);

disp('SPM segment finished');
BaselineNoBadSliFN2=[WorkingP 'BaselineNoBadSli2.nii'];

% In unix, run the system cp command with no "-p" because it gives an
% error when the destination is in another computer so source and dest
% files have different owner
if (filesep == '/') % Unix
    system(['cp -f ' BaselineNoBadSliFN ' ' BaselineNoBadSliFN2]);
else  % Windows
    copyfile(BaselineNoBadSliFN,BaselineNoBadSliFN2,'f');    
end

% DCEMeanSegP2=SPM_Segment(BaselineNoBadSliFN2,Force,[],FMaskFN);
DCEMeanSegP2=DCEMeanSegQB1;

BaseSeg3D(:,:,:,1)=loadniidata([DCEMeanSegP2 'c1ForSeg.nii'])/256;
BaseSeg3D(:,:,:,2)=loadniidata([DCEMeanSegP2 'c2ForSeg.nii'])/256;
BaseSeg3D(:,:,:,3)=loadniidata([DCEMeanSegP2 'c3ForSeg.nii'])/256;
BaseSeg3D(:,:,:,4)=loadniidata(FMaskFN);
[Tmp, BaseSeg3DAll]=max(BaseSeg3D(:,:,:,1:3),[],4);
BaseSeg3DAll(~BaseSeg3D(:,:,:,4))=0;
% BaseCleaned=loadniidata([DCEMeanSegP2 'mForSeg.nii']);
BaseCleaned=MeanVol;

se=strel('disk',4,8);
EBrainMask=imerode(FBrainMask,se);
EEBrainMask=imerode(EBrainMask,se);
EEEBrainMask=imerode(EEBrainMask,se);

[MaxVal3D MaxEnhancementTime]=max(DCE4D,[],4);
MaxRatioToBaseline=MaxVal3D./Baseline;
EnhancementMsk=MaxRatioToBaseline>MinEnhancementR & MaxEnhancementTime>(BolusStart-1); % & MRIdx<(BolusStart+3);
TooEnhancedForNAWM=MaxRatioToBaseline>1.5;

CSFMask=BaseSeg3DAll==3 & BaseSeg3D(:,:,:,3)>Options.ThreshForRefMasks & EEEBrainMask==1;
WMMask=BaseSeg3DAll==2 & BaseSeg3D(:,:,:,2)>Options.ThreshForRefMasks & EEEBrainMask==1 & ~TooEnhancedForNAWM;

if(sumn(WMMask)<100)
    WMMask=BaseSeg3DAll==2 & BaseSeg3D(:,:,:,2)>Options.ThreshForRefMasks*0.9 & EEEBrainMask==1 & ~TooEnhancedForNAWM;
end
BaseSeg3DAllx=BaseSeg3DAll;
BaseSeg3DAllx(CSFMask)=4;
BaseSeg3DAllx(WMMask)=5;

Tmp=max(BrainMask,[],3);
F=find(max(Tmp,[],2));
GoodRows=F(1):F(end);
F=find(max(Tmp,[],1));
GoodCols=F(1):F(end);
GoodSlices=setdiff(1:SDCE(3),BadSlicesF2);
%%
%     MaxV=max(BaseCleaned(:))*0.1;
MaxV=median(BaseCleaned(isfinite(BaseCleaned) & BaseCleaned>100))*2;
[q MaxV]=FindDR(BaseCleaned(BrainMask));
for i=1:numel(GoodSlices)
    CurSli=GoodSlices(i);
    I=squeeze(BaseCleaned(:,:,CurSli));
    %         [Tmp MaxV]=FindDR(I);
    ClrM=[1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1];
    IRGB=repmat(min(1,I/MaxV),[1 1 3]);
    for tt=1:3
        BW2 = bwmorph(squeeze(BaseSeg3DAll(:,:,CurSli)==tt),'remove');
        for kk=1:3
            TmpI=IRGB(:,:,kk);
            TmpI(BW2)=ClrM(tt,kk);
            IRGB(:,:,kk)=TmpI;
        end
    end
    for tt=4:5
        BW2 = bwmorph(squeeze(BaseSeg3DAllx(:,:,CurSli)==tt),'remove');
        for kk=1:3
            TmpI=IRGB(:,:,kk);
            TmpI(BW2)=ClrM(tt,kk);
            IRGB(:,:,kk)=TmpI;
        end
    end
    IRGB3(:,:,:,i)=IRGB;
end

figure(9899);clf;
montage(mritransform(IRGB3(GoodRows,GoodCols,:,:)))
title(num2str(GoodSlices));
saveas(9899,[WorkingP 'BaseSeg'  '.png']);
saveas(9899,[WorkingP 'BaseSeg'  '.fig']);
%%
close(9899);
AddToLog(WorkingP,'ya_2d',['Img segmentation. Red - GM, Green - WM, Blue - CSF, Magenta - WM for reference, Yellow - CSF for reference.'],['BaseSeg'  '.png']);
%%
Raw2Nii(CSFMask,[WorkingP 'RefAuto_Base' '_CSF_2430.nii'],'float32', MeanFN);
Raw2Nii(WMMask,[WorkingP 'RefAuto_Base' '_WM_830.nii'],'float32', MeanFN);
%% Step 7 - Enhancement mask
% disp('SPM segmentation');
NoEnhancementVol=Baseline;
CNoEnhancementVol=NoEnhancementVol;
CNoEnhancementVol(~BrainMask)=0;
CNoEnhancementFN=[WorkingP 'CleanedNoEnhancementNo' GroupToStr(num2strC(BadSlicesF2)) '.nii'];
Raw2Nii(CNoEnhancementVol,CNoEnhancementFN,'float32', MeanFN);
% DCEMeanSegP=SPM_Segment(CNoEnhancementFN,Force,[],false,BrainMskFN);
% %%
% DCEMeanSegP=SPM_Segment(MeanFNGoodSlices,Force,[],false,BrainMskFN);
%% Step 8 - Final mask out of all masks + save data

Msk=EnhancementMsk & MskSlices & MskMinSignal & BrainMask;
% Get the most interesting slice
[QQ, RepSli]=max(squeeze(gsum(Msk,1:2)));
% DCE2D=Reshape4d22d(DCE4D,Msk);
% Baseline2D=Baseline(Msk);
% RDCE2D=repMulti(sDCE2D,1./Baseline2D);

%% Bad time points check

% Start point to check will be 10 volumes after the bolus start volume
StartPointForBad=BolusStart+10;
% Get the median plot from start point to end
median_to_end=(MedTC(StartPointForBad:end));
% Smooth the data
[dataout lowerLimit upperLimit smooth_only] = lowess([1:numel(median_to_end);median_to_end]',0.1,0);
Smoothed=dataout(:,3)';
% Error pf median to smoothed
DiffOrigSmoothed=(median_to_end-Smoothed);
% Get a normalized error vector
RobustStd=median(abs(DiffOrigSmoothed))/0.67;
% NDiffOrigSmoothed=(DiffOrigSmoothed-mean(DiffOrigSmoothed))/std(DiffOrigSmoothed);
NDiffOrigSmoothed=(DiffOrigSmoothed-mean(DiffOrigSmoothed))/RobustStd;
% Points with higher value than threshold will be considered bad time points
ThreshForBadTimePoints=2.5;
% Get the volume index of the bad points
BadTimePoints=find(abs(NDiffOrigSmoothed)>ThreshForBadTimePoints)+StartPointForBad-1;
% Create graphs displaying the finding bad points process
figure(121);clf;subplot(2,1,1);
plot(1:nVols,MedTC,'r.',StartPointForBad:nVols,Smoothed,'b-',BadTimePoints,MedTC(BadTimePoints),'ko');
title('Original and smoothed median time course');
subplot(2,1,2);
plot(1:nVols,[zeros(1,StartPointForBad-1) NDiffOrigSmoothed],'g-',1:nVols,(1:nVols)*0+ThreshForBadTimePoints,'r-',1:nVols,(1:nVols)*0-ThreshForBadTimePoints,'r-',BadTimePoints,NDiffOrigSmoothed(BadTimePoints-StartPointForBad+1),'ko');
title(['Normalized diff, bad: ' num2str(BadTimePoints)]);
saveas(121,[WorkingP 'BadTimePoints.fig']);
saveas(121,[WorkingP 'BadTimePoints.png']);
close(121);
if(~isempty(BadTimePoints))
    Txt=['Bad time points detected: ' num2str(BadTimePoints) '.'];
else
    Txt='No bad time points detected.';
end
AddToLog(WorkingP,'a_2ff',Txt,[],1);
AddToLog(WorkingP,'ya_2ff',Txt,'BadTimePoints.png');
%% show first vol
Tmp=max(Msk,[],3);
F=find(max(Tmp,[],2));
GoodRows=F(1):F(end);
F=find(max(Tmp,[],1));
GoodCols=F(1):F(end);
%Tmp=loadniidata([WorkingP '\DCEMainCoreged\Coreged_vol_0001.nii']);
% Tmp=loadniidata([WorkingP filesep 'DCEMainCoreged' filesep 'Coreged_vol_0001.nii']);
Tmp=loadniidata(DCEFNs{1});
Tmp=Tmp(GoodRows,GoodCols,:);
[QQ, MX]=FindDR(Tmp(:));
Tmp=permute(Tmp,[1 2 4 3])./MX;
Tmp=repmat(Tmp,[1 1 3 1]);
figure(1111);clf;
montage(mritransform(Tmp));
title('First DCE vol.');
saveas(1111,[WorkingP 'FirstVol.png']);
saveas(1111,[WorkingP 'FirstVol.fig']);
close(1111);
clear Tmp
AddToLog(WorkingP,'xa_2firstvol','First volume of DCE.','FirstVol.png');
%%
DCE4DX=max(0,DCE4D-repmat(Baseline,[1 1 1 SDCE(4)]));
DCE4DX(:,:,BadSlicesF,:)=NaN;
[MaxVal, PeakTime]=max(DCE4DX,[],4);
PeakTime(~FBrainMask)=NaN;
GoodPeak=PeakTime>(BolusStart-3) & PeakTime<(BolusStart+3);
[I J K]=ind2sub(size(GoodPeak),find(GoodPeak));
NearVals=[DCE4DX(sub2ind(size(DCE4D),I,J,K,PeakTime(GoodPeak)-1)) DCE4D(sub2ind(size(DCE4D),I,J,K,PeakTime(GoodPeak)+1))];
MaxVal1D=MaxVal(GoodPeak);
PeakTime1D=PeakTime(GoodPeak);
[MaxNearVal, WhichDirection]=max(NearVals,[],2);
RelativeDist=MaxNearVal./(MaxVal1D+MaxNearVal);
MoreExactPeakTime1D=PeakTime(GoodPeak)+((WhichDirection-1)*2-1).*RelativeDist;
MoreExactPeakTime=PeakTime;
MoreExactPeakTime(GoodPeak)=MoreExactPeakTime1D;
% MoreExactPeakTime(~GoodPeak)=NaN;
clear DCE4DX;
%%
save(PrepareFN);
clear DCE4D DCE2D Baseline2D
saveas(78362,[WorkingP 'SlicesIntensityAndBolusTime.png']);
saveas(78362,[WorkingP 'SlicesIntensityAndBolusTime.fig']);
close(78362)
AddToLog(WorkingP,'a_2xx',['Bolus start index estimation:' num2str(BolusStart) '.']);
AddToLog(WorkingP,'ya_2xx',['Bolus start index estimation:' num2str(BolusStart) '.'],'SlicesIntensityAndBolusTime.png');
disp(['DCE_Prepeare4D, end ' datestr(now)]);
% %%
% Tmp=Baseline+BrainMaskA*500;
% 
% [QQ, MX]=FindDR(Tmp(:));
% Tmp=permute(Tmp,[1 2 4 3])./MX;
% Tmp=repmat(Tmp,[1 1 3 1]);
% figure(1111);clf;
% montage(mritransform(Tmp));
% 
% % figure;imagesc(mritransform(Tmp(:,:,MidSli)))