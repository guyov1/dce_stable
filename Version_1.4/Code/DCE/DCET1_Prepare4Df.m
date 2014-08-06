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

% ASK GILAD - what is the purpose of the following?
% ASNWER - the time between each slice is dependent on the resolution, repetition time and acceleration ingredients.
TimeBetweenDCEVolsPerSlice=double(min(CurFullInfo.AcquisitionMatrix(CurFullInfo.AcquisitionMatrix>0))*CurFullInfo.RepetitionTime)/(1000*max(CurFullInfo.EchoTrainLength,1))*CurFullInfo.PercentPhaseFieldOfView/100.;
%TimeBetweenDCEVolsPerSlice=double(min(CurFullInfo.AcquisitionMatrix(CurFullInfo.AcquisitionMatrix>0))*CurFullInfo.RepetitionTime)/(1000*max(CurFullInfo.EchoTrainLength,1));

%% Check if stage was done already / "force" was used.

% Path for the output data calculated in this stage
PrepareFN=[WorkingP 'AfterPrepare4D.mat'];

% If the data already exists, don't calculate again
if(exist(PrepareFN,'file') && ~Force)
    D=dir(PrepareFN);
%     if(D.datenum>DDate)
%         loadBut(PrepareFN,{CurVars.name});
        disp('DCE_Prepeare4D already computed');
        return;
%     end
end

% If the file exists and we got here, we used Force -> Delete anyway.
if(exist(PrepareFN,'file'))
    delete(PrepareFN);
end

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
DCE4D=loadCoregedNiftis(DCEFNs);
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

%% Step 3 - Create mutual information curve to see when the registration did not succeed

disp('Calculate mutual information');
% Mutual information curve
% ASK GILAD - explain how he calculated the MI and what it means exactly.
MIC=CalcMICurve(DCEFNs,MeanFN);

%%  Step 4 - Show movements -> 	Create a montage of all middle slices during time.

disp('Disply coreged');
% if(false)
    % Middle Slice
    MidSli=floor(SDCE(3)/2);
    
    % Dynamic Range?
    % ASK GILAD - why does he divide by 3?
    DR=[min(DCE4D(:)) max(DCE4D(:)/3)];
        figure(987);clf;
        
    % Cretae a montage image of all the middle slices during time.
    montage(mritransformNoSqueeze((DCE4D(:,:,MidSli,:))-DR(1))./(DR(2)-DR(1)));
    gprint(987,[WorkingP 'CoregedMidSlice.jpg']);
    close(987);
% end

%% Step 5 -  Masks creation

% TimeBetweenDCEVols
TimeBetweenDCEVols=TimeBetweenDCEVolsPerSlice*SDCE(3);
TimeBetweenDCEVolsMin=TimeBetweenDCEVols/60;

% Frequency of images over time
% The following variable is used later in DCET1_CTCf.m
Min2Timestep=1/TimeBetweenDCEVolsMin;

%% Min signal mask

%  Mask all time volumes to see whether in the entire test, the signal was stronger than the minimal.
MskMinSignal=all(DCE4D>MinSignal,4);

%% Masking by FSL (previously used SPM)
disp('Brain extraction');

% Running FSL's bet function to extract the brain out of the image (mask)
[BetStrippedMeanFN, BetMaskFN]=mcbet2(MeanFN,Force);
BrainMask=loadniidata(BetMaskFN)>0;

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
BadSlices(~NaNSlices)=abs(MedSlice(~NaNSlices)-mean(MedSlice(MidSli-1:MidSli+1)))>50;
BadSlicesF=find(BadSlices);
disp(['Ignoring slices ' num2str(BadSlicesF')]);
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

%% Step 7 - Enhancement mask
% disp('SPM segmentation');

NoEnhancementVol=Baseline;
CNoEnhancementVol=NoEnhancementVol;
CNoEnhancementVol(~BrainMask)=0;
CNoEnhancementFN=[WorkingP 'CleanedNoEnhancementNo' GroupToStr(num2strC(BadSlicesF2)) '.nii'];
Raw2Nii(CNoEnhancementVol,CNoEnhancementFN,'float32', MeanFN);
% DCEMeanSegP=SPM_Segment(CNoEnhancementFN,Force,[],false,BrainMskFN);
% %%
% % DCEMeanSegP=SPM_Segment(MeanFNGoodSlices,Force,[],false,BrainMskFN);
% C1=loadniidata([DCEMeanSegP 'c1ForSeg.nii'])/256;
% C2=loadniidata([DCEMeanSegP 'c2ForSeg.nii'])/256;
% C3=loadniidata([DCEMeanSegP 'c3ForSeg.nii'])/256;
% MinSPMBrainValue=0.3;
% % BrainMask=(C1+C2)>MinSPMBrainValue;
% % 
% % BrainMaskA=(C1+C2+C3)>MinSPMBrainValue;
% BrainMaskA=bwfillHoles3Dby2D(BrainMask);
% disp('SPM segment finished');


[MaxVal3D MaxEnhancementTime]=max(DCE4D,[],4);
MaxRatioToBaseline=MaxVal3D./Baseline;
EnhancementMsk=MaxRatioToBaseline>MinEnhancementR & MaxEnhancementTime>(BolusStart-1); % & MRIdx<(BolusStart+3);

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
NDiffOrigSmoothed=(DiffOrigSmoothed-mean(DiffOrigSmoothed))/std(DiffOrigSmoothed);
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
saveas(121,[WorkingP 'BadTimePoints.jpg']);
close(121);

%%
% Output to text
fid=fopen([WorkingP 'Preprocessdata.csv'],'w');
fprintf(fid,'TimeBetweenDCEVols,%f\n',TimeBetweenDCEVols);
fprintf(fid,'BadTimePoints,',TimeBetweenDCEVols);
for i=1:numel(BadTimePoints)
    fprintf(fid,'%d ',BadTimePoints(i));
end
fprintf(fid,'\n');
fclose(fid);
save(PrepareFN);
clear DCE4D DCE2D Baseline2D
gprint(78362,[WorkingP 'SlicesIntensityAndBolusTime.jpg']);
close(78362)
disp(['DCE_Prepeare4D, end ' datestr(now)]);