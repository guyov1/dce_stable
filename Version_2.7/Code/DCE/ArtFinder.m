AddToLog(WorkingP,'cx_00','\\subsection*{Finding the Arteries}');
if(~exist('DebugSuffix','var'))
    DebugSuffix='';
end

%%
% load(RelaxFN,'GoodSlices','GoodRows','GoodCols');
Tmp=max(FBrainMask,[],3);
F=find(max(Tmp,[],2));
GoodRows=F(1):F(end);
F=find(max(Tmp,[],1));
GoodCols=F(1):F(end);
load(PrepareFN,'BadSlicesF2');
GoodSlices=setdiff(1:size(FBrainMask,3),BadSlicesF2);


Msk3D=Idx3D>0;
MskE=Msk3D;
MskCTC=Msk3D;
FBrainMask=bwfillHoles3Dby2D(Msk3D);

% Erodes each of the slices
for i=1:size(Msk3D,3)
    MskE(:,:,i)=imerode(squeeze(FBrainMask(:,:,i)),strel('disk',16));
end
%% Filters
MskNotInEdge=MskE & Msk3D;
[DataOutA, ~, OtherMasks, Msks]=FilterCTCToFindArt(CTC2D(:,1:NumVols),Msk3D,NumVols,TimeBetweenDCEVolsFinal,BolusStart,Options);
for i=1:numel(Msks)
    Msks{i}=Msks{i} & MskNotInEdge;
    MsksN(i)= sumn(Msks{i});
end
AllMsks=cat(4,Msks{:});
SumMsk=sum(AllMsks,4);
NinNsks=histc(SumMsk(:),0:numel(Msks));
TrgN=300;
if(NinNsks(end)>TrgN)
    MskF=SumMsk==numel(Msks);
else
    MskAI=numel(Msks)-find(cumsum(flipud(NinNsks))>TrgN,1)+1;
    MskA=SumMsk>MskAI;
    MskB=(SumMsk==MskAI);
    TmpF=find(MskB);
    Needed=TrgN-sumn(MskA);
    RIdxs=randsample(numel(TmpF),Needed);
    Tmp=TmpF*0>1;
    Tmp(RIdxs)=true;
    MskB(MskB)=Tmp;
    MskF=MskA | MskB;
end
% MaskWithEnoughI=find(MsksN>10,1,'last');
% MskF=Msks{MaskWithEnoughI} & MskNotInEdge;
% MskF=OtherMasks & MskNotInEdge;
InIdxs=getIndicesOfMskInsideMsk(MskF,Msk3D);
DataOut=CTC2D(InIdxs,:);
Raw2Nii(MskF,[WorkingP 'ArtMsk' DebugSuffix '.nii'],'float32',MeanFN);
%%
ManualArtMaskFN=[WorkingP 'ManualArtMask' DebugSuffix '.nii'];
if(strcmp(ROIStr,'Full'))
    if(exist(ManualArtMaskFN,'file'))
        ManualArtMask=loadniidata(ManualArtMaskFN);
        MskF=MskF & ManualArtMask;
        InIdxs=getIndicesOfMskInsideMsk(MskF,Msk3D);
        DataOut=CTC2D(InIdxs,:);
    end
else
    tmp=getKthElement(getComputerParams('tpm'),1);
    ROIsPath=[tmp(1:end-23) 'Code' filesep 'Utils' filesep 'SPM_precofigures' filesep];
    SrcROI=[ROIsPath ROIStr '.nii'];
    TrgFN=NormalizeWrite(SrcROI,[WorkingP 'Seg_QuickB1Corrected_DCEMean' filesep 'ForSeg_seg_inv_sn.mat'],true,WorkingP);
    rTrgFN=CoregWrite(TrgFN,getTMatFromNii(TrgFN),true,WorkingP,false,MeanFN);
    %     mricronx({MeanFN TrgFN})
    ManualArtMask=loadniidata(rTrgFN);
    MskF=MskF & ManualArtMask;
    InIdxs=getIndicesOfMskInsideMsk(MskF,Msk3D);
    DataOut=CTC2D(InIdxs,:);
end

%%  Get the representing voxels data and noise
% [CVI, BinCVI, Bin2CVI DataToFitA]=ChooseRepVoxelsForAIFFind(Msk3D,MskE,CTC2D(:,1:NumVols),BolusStart,Options.Rep_MaxAroundBolus,Options.Rep_RatioToEnd,Options.Rep_nPerSet);
[CVI, BinCVI, Bin2CVI DataToFitA]=ChooseRepVoxelsForAIFFind(MskF,MskE,DataOut,BolusStart,Options.Rep_MaxAroundBolus,Options.Rep_RatioToEnd,Options.Rep_nPerSet);
CVI=InIdxs(CVI);
% Options=struct('SubSecRes',{[]},'MaxTDif',{[]});
%Get the relative noise (percentage out of max signal value ) of each representing voxel( rmadCTC2D is the noise of each one)
% DataNoise=rmadCTC2D(CVI);
DataNoise=EstimateNoise(DataToFitA);
% The data ( c(t) ) of the representing voxels
DataToFit=CTC2D(CVI,:);
DataToFit=DataToFitA;
% Plot2DDataOnSubfigures(11,DataToFit);

DataToFit2=DataToFit;
DataNoise=EstimateNoise(DataToFit2);
CVI2=CVI;
BinCVI2=BinCVI;
Bin2CVI2=Bin2CVI;
MskCTC2=MskCTC;
CTC2D2=CTC2D;
% Plot2DDataOnSubfigures(12,DataToFit2)
%%
ShowData=false;
if(ShowData)
    gfig(15);
    subplot(1,2,1);
    plot(DataToFit2');
    subplot(1,2,2);
    plot(NormalizeByRows(DataToFit2)');
    hold on;
    plot(mean(NormalizeByRows(DataToFit2))','k','LineWidth',2)
end
%%
GoodTs=Sts;
GoodTIdxs=1:numel(Sts);
if (Num_Of_T1_Maps>1)
    NumBefore=size(Additional_before_main,2);
    AfterIdxs=find(Additional_T1_Maps_Time_Diff_Sets>0)';
%     GoodTIdxs=[(NumBefore+1):(numel(ts)+NumBefore) numel(ts)+SingleAngleIdxs-1];
%     GoodTs=[Sts(1:numel(ts)) Additional_T1_Maps_Time_Diff_Sets(SingleAngleIdxs)'*86400/60];
    
    GoodTIdxs=[(NumBefore+1):(numel(ts)+NumBefore) numel(ts)+AfterIdxs-1];
    GoodTs=[Sts(1:numel(ts)) Additional_T1_Maps_Time_Diff_Sets(AfterIdxs)'*86400/60];
    
    GoodTIdxs=1:numel(SampleTsNoBefore);
    GoodTs=SampleTsNoBefore/60;
    % figure;plot(GoodTs,DataToFit(1:10:end,GoodTIdxs)','.-');
end

RepVoxFN=[WorkingP 'RepVox' USStr DebugSuffix '.mat'];

Mx=max(CTC2D,[],2);
SMx=sort(Mx);
MaxAmp=SMx(numel(Mx)-10);
% Show the Chosen
%     Tmp=ShowImageWithPoints(1122,CT1,MskCTC,CVI,GoodRows,GoodCols,GoodSlices);
Tmp=ShowImageWithPoints(1122,CT1,MskCTC2,CVI2,GoodRows,GoodCols,GoodSlices);
Raw2Nii(Tmp*1,[WorkingP 'ChosenVoxelsForAIFFinding' DebugSuffix '.nii'],'float32', MeanFN);
RelaxFN=[WorkingP 'Relax.mat'];
%%
title('Chosen for AIF estimation');
saveas(1122,[WorkingP 'ChosenForAIFEstimation' DebugSuffix '.png']);
saveas(1122,[WorkingP 'ChosenForAIFEstimation' DebugSuffix '.fig']);
close(1122);
AddToLog(WorkingP,'ycx_at1','Chosen for AIF estimation','ChosenForAIFEstimation.png');
%% 
AIFFinderFN=[WorkingP 'AIFFindData' regexprep(USStr,'\.','_') DebugSuffix '.mat'];

clear CTC2DBigGood B1Maps OAdditional_T1_Maps AdditionalT1_Inverse_2D FinalT1s adCTC2D B1 BaselineB1Fldx BaselineFA3Dx

% Save the stage .mat file
save(RepVoxFN,'DataToFit2','CVI2','MskCTC2','CTC2D2');