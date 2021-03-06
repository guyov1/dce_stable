function DCET1_Prepare4Df(CurMainDCEInfo,WorkingP,Force, Options)

CurMainDCEInfo.Path=strrep(CurMainDCEInfo.Path,'//','/');
DCEMainP=CurMainDCEInfo.Path;
CurFullInfo=dicominfo(CurMainDCEInfo.Filename);
TimeBetweenDCEVolsPerSlice=double(min(CurFullInfo.AcquisitionMatrix(CurFullInfo.AcquisitionMatrix>0))*CurFullInfo.RepetitionTime)/(1000*max(CurFullInfo.EchoTrainLength,1));

%% Check in longitudinal
PrepareFN=[WorkingP 'AfterPrepare4D.mat'];
if(exist(PrepareFN,'file') && ~Force)
    D=dir(PrepareFN);
%     if(D.datenum>DDate)
%         loadBut(PrepareFN,{CurVars.name});
        disp('DCE_Prepeare4D already computed');
        return;
%     end
end
if(exist(PrepareFN,'file'))
    delete(PrepareFN);
end
a=rmdir([WorkingP 'Seg_DCEMean'],'s');
a=rmdir([WorkingP 'Seg_DCEMeanL'],'s');
a=rmdir([WorkingP 'DCEMainCoregedLongitudinal'],'s');

disp(['DCE_Prepeare4D, start ' datestr(now)]);
MinSignal=50;
MinEnhancementR=1.02;
% MinSPMBrainValue=0.3;
%% Step 1 - Clean and realign DCE images
disp('Dicom to Nifti');
DCEMNiiOrigP=[WorkingP 'DCEMainNii' filesep];
mkdir(DCEMNiiOrigP);
DDCE=dir([DCEMNiiOrigP '*.nii']);
if(isempty(DDCE) || Force)
    gDicom2Nifti(DCEMainP,[DCEMNiiOrigP 'vol.nii']);
    
    DDCE=dir([DCEMNiiOrigP '*.nii']);
    if(isempty(DDCE))
        error('DCET1_Prepare4Df - Dicom2Nifi didnt work');
    end
end
% Clean - not needed?
% Realign and extract mean
disp('Coregistration');
DCE_Coreg;
DCEMNiiP=DCECoregP;
% if(CoregLongPerformed)
%     DCE_CoregLong;
% end
%% Load all and create 4D and mean files
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

%%
disp('Calculate mutual information');
% Mutual information curve
MIC=CalcMICurve(DCEFNs,MeanFN);
%% Show movements
disp('Disply coreged');
% if(false)
    MidSli=floor(SDCE(3)/2);
    DR=[min(DCE4D(:)) max(DCE4D(:)/3)];
    figure(987);clf;
    montage(mritransformNoSqueeze((DCE4D(:,:,MidSli,:))-DR(1))./(DR(2)-DR(1)));
    gprint(987,[WorkingP 'CoregedMidSlice.jpg']);
    close(987);
% end
%% TimeBetweenDCEVols
TimeBetweenDCEVols=TimeBetweenDCEVolsPerSlice*SDCE(3);
TimeBetweenDCEVolsMin=TimeBetweenDCEVols/60;
Min2Timestep=1/TimeBetweenDCEVolsMin;
%%%%%% Masks
%% Min signal
MskMinSignal=all(DCE4D>MinSignal,4);
%% Masking by SPM
disp('Brain extraction');

BadSlicesA=[1 size(MeanVol,3)];
% MeanVolGoodSlices=MeanVol;
% MeanVolGoodSlices(:,:,BadSlicesA)=0;
% MeanFNGoodSlices=[WorkingP 'MeanFNNo' GroupToStr(num2strC(BadSlicesA)) '.nii'];
% Raw2Nii(MeanVolGoodSlices,MeanFNGoodSlices,'float32', MeanFN);
% DCEMeanSegP=SPM_Segment(MeanFNGoodSlices,Force,[],false);
% 
% % DCEMeanSegP=SPM_Segment(MeanFN,false,[],false);
% C1=loadniidata([DCEMeanSegP 'c1ForSeg.nii'])/256;
% C2=loadniidata([DCEMeanSegP 'c2ForSeg.nii'])/256;
% BrainMask=(C1+C2)>MinSPMBrainValue;
% SPMC1C2MskFN=[WorkingP 'SPMC1C2Msk.nii'];
% Raw2Nii(BrainMask,SPMC1C2MskFN,'float32',MeanFN);

[BetStrippedMeanFN, BetMaskFN]=mcbet2(MeanFN,Force);
BrainMask=loadniidata(BetMaskFN)>0;
BrainMask(:,:,BadSlicesA)=false;
% BrainMask=gErode(BrainMaskA,[5 5 1]);
% %% Coregister T1Map to mean DCEMain
% T1FN=[WorkingP 'T1.nii'];
% CrgMatFN=CoregEstimate(T1FN,MeanFN);
% CT1FN=CoregWrite(T1FN,CrgMatFN);
%% Selecting slices
disp('Find problematic slices');

MedSlice=zeros(SDCE(3),1);
for i=1:SDCE(3)
    SliceMsk=zeros(SDCE(1),SDCE(2),SDCE(3));
    SliceMsk(:,:,i)=BrainMask(:,:,i);
    Tmp=Reshape4d22d(DCE4D(:,:,:,1),SliceMsk);
    MedSlice(i,:)=median(Tmp,1);
end
% NMedSlice=repMulti(repPlus(MedSlice,-mean(MedSlice,1)),1./std(MedSlice,0,1));
% BadSlices=any(abs(NMedSlice)>2,2);
NaNSlices=isnan(MedSlice);
disp('Looking for weird slices');
if(numel(MedSlice(~NaNSlices))<2)
    error('Only one slice!');
end
% [QQ, optmixture] = GaussianMixture(MedSlice(~NaNSlices), 3, 0,false);
% [QQ, BigGroupI]=max([optmixture.cluster.pb]);
% OthergroupsI=setdiff(1:optmixture.K,BigGroupI);
BadSlices=NaNSlices;
% BadSlices(~NaNSlices)=sum(optmixture.pnk(:,OthergroupsI),2)>0.5;
BadSlices(~NaNSlices)=abs(MedSlice(~NaNSlices)-mean(MedSlice(MidSli-1:MidSli+1)))>50;
BadSlicesF=find(BadSlices);
disp(['Ignoring slices ' num2str(BadSlicesF')]);
figure(78362);clf;subplot(1,2,1);
plot(1:SDCE(3),MedSlice,'b',BadSlicesF,MedSlice(BadSlicesF),'ro');
title('Median of slices and ignored ones');
xlabel('Slices');
ylabel('Median');
MskSlices=ones(size(MskMinSignal))>0;
MskSlices(:,:,BadSlicesF)=false;
%%
BadSlicesF2=union(BadSlicesF,[1 size(MeanVol,3)]);
BrainMskFN=[WorkingP 'BrainMask.nii'];
BrainMask(:,:,BadSlicesF2)=0;
Raw2Nii(BrainMask,BrainMskFN,'float32',MeanFN);
MeanVolGoodSlices=MeanVol;
MeanVolGoodSlices(~BrainMask)=0;
MeanVolGoodSlices(:,:,BadSlicesF2)=0;
MeanFNGoodSlices=[WorkingP 'MeanFNNo' GroupToStr(num2strC(BadSlicesF2)) '.nii'];
Raw2Nii(MeanVolGoodSlices,MeanFNGoodSlices,'float32', MeanFN);
%% Find bolus start, compute baseline
disp('Rough estimation of bolus time');

DCE2D=Reshape4d22d(DCE4D,MskMinSignal);
MedTC=median(DCE2D,1);

% [~, optmixture] = GaussianMixture(MedTC', 3, 0,false);
% [Z Grouped]=max(optmixture.pnk,[],2);
% if(max(Grouped)==1)
%     Smoothed=conv2(MedTC,ones(1,ceil(nVols/10)),'same');
%     [a BolusStart]=max(Smoothed);
% else
%     Grouped(end)=Grouped(1)+1;
%     BolusStart=find(Grouped~=Grouped(1),1);
% end
% Second method
TwoMinTimePoint=floor(2/TimeBetweenDCEVolsMin);
Ps=zeros(1,numel(MedTC))+2;
for i=3:min(TwoMinTimePoint,numel(MedTC)-2)
    [h Ps(i)]=ttest2(MedTC(1:i),MedTC((i+1):end),[],[],'unequal');
end
mLPs=-log(Ps);
% figure;plot(1:numel(MedTC),MedTC,'b',1:numel(MedTC),mLPs.*(max(MedTC)-min(MedTC))./(max(mLPs)-min(mLPs))+min(MedTC),'r')
[Tmp, BolusStart]=max(mLPs);
BolusStart=BolusStart+1;
BolusStartMin=(BolusStart-1)*TimeBetweenDCEVolsMin;
% BolusStart=find(MedTC>MedTC(1)+20,1);
Baseline=mean(DCE4D(:,:,:,1:(BolusStart-2)),4);
BaselineFN=[WorkingP 'Baseline.nii'];
Raw2Nii(Baseline,BaselineFN,'float32', MeanFN);
figure(78362);subplot(1,2,2);
plot(MedTC); hold on;plot([BolusStart BolusStart],[min(MedTC) max(MedTC)],'r');
title('Bolus start approximation');
%%
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
%% Enhancement mask
[MaxVal3D MaxEnhancementTime]=max(DCE4D,[],4);
MaxRatioToBaseline=MaxVal3D./Baseline;
EnhancementMsk=MaxRatioToBaseline>MinEnhancementR & MaxEnhancementTime>(BolusStart-1); % & MRIdx<(BolusStart+3);
%% Final mask
Msk=EnhancementMsk & MskSlices & MskMinSignal & BrainMask;
[QQ, RepSli]=max(squeeze(gsum(Msk,1:2)));
% DCE2D=Reshape4d22d(DCE4D,Msk);
% Baseline2D=Baseline(Msk);
% RDCE2D=repMulti(DCE2D,1./Baseline2D);
%%
save(PrepareFN);
clear DCE4D DCE2D Baseline2D
gprint(78362,[WorkingP 'SlicesIntensityAndBolusTime.jpg']);
close(78362)
disp(['DCE_Prepeare4D, end ' datestr(now)]);