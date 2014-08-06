function DCET1_RelaxForSubjectf(DCEInfos,WorkingP,DoN3,DoGlobal,DoDeviations,CalcForce,WhichMean,Options)

DCERelaxP=[WorkingP 'Relaxometry' filesep];
RWorkingP=DCERelaxP;
if(exist([RWorkingP 'T13DNFA_N3k.nii'],'file') && ~CalcForce)
    disp('DCET1_RelaxForSubjectf Already computed');
    return;
end
disp('Starting DCET1_RelaxForSubjectf');

TRsF=[DCEInfos.RepetitionTime];
FAsF=[DCEInfos.FlipAngle];
%% Base
CMeanFN=[WorkingP 'CopyDCEMean.nii'];
CurOtherFAsInfos=DCEInfos;
nOtherFAs=numel(CurOtherFAsInfos);
OtherFACrgFNs=cell(1,nOtherFAs);
OtherFABaseFNs=cell(1,nOtherFAs);
OtherFAMatFNs=cell(1,nOtherFAs);
mkdir(DCERelaxP);
%% Dicom2Nifti
for i=1:nOtherFAs
    OtherFABaseFNs{i}=[DCERelaxP 'OrigFA_' sprintf('%02d',CurOtherFAsInfos(i).FlipAngle) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '.nii'];
    Tmp=gDicom2Nifti(CurOtherFAsInfos(i).Path,OtherFABaseFNs{i});
    MainSubFN=[DCERelaxP 'OrigFA_' sprintf('%02d',CurOtherFAsInfos(i).FlipAngle) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '_01.nii'];
    if(Tmp>1 && Tmp<10 && exist(MainSubFN,'file'))
        movefile(MainSubFN,OtherFABaseFNs{i});
        delete([DCERelaxP 'OrigFA_' sprintf('%02d',CurOtherFAsInfos(i).FlipAngle) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '_*.nii']);
    end
    if(~exist(OtherFABaseFNs{i},'file'))
        error('DCET1_RelaxForSubjectf Dcm2Nii NotCreated');
    end
end
%% Coreg dest -
MedFA=median(FAsF);
MedFAI=find(FAsF==MedFA,1);
if(strcmp(WhichMean,'Mean 4D'))
    MeanFN=[WorkingP 'DCEMean.nii'];
else
    if(strcmp(WhichMean,'Median angle'))
        MeanFN=OtherFABaseFNs{MedFAI};
    else
        MeanFN=WhichMean;
    end
end
delete(CMeanFN);
try
    copyfile(MeanFN,CMeanFN,'f');
catch MErr
end
if(~exist(CMeanFN,'file'))
    error('DCET1_RelaxForSubjectf: Problem copying CMeanFN');
end
%% Coreg
for i=1:nOtherFAs
    OtherFAMatFNs{i}=CoregEstimate(OtherFABaseFNs{i},CMeanFN,CalcForce);
    OtherFACrgFNs{i}=CoregWrite(OtherFABaseFNs{i},OtherFAMatFNs{i},CalcForce,DCERelaxP,false,CMeanFN);
end
delete(CMeanFN);
disp('Relaxometry coregistration finished');
%% Full brain mask for relaxometry
PrepareFN=[WorkingP 'AfterPrepare4D.mat'];

load(PrepareFN,'BrainMask','C2','RepSli','MeanVol');
FBrainMask=bwfillHoles3Dby2D(BrainMask);
A=load_untouch_nii(OtherFACrgFNs{1});
A.img=int16(FBrainMask);
FMaskFN=[DCERelaxP 'FBrainMsk.nii'];
save_untouch_nii(A,FMaskFN);
disp('FBrainMsk finished');
%% Segmentation for masks
se=strel('disk',4,8);
EBrainMask=imerode(FBrainMask,se);
if(DoGlobal)
    CCSplenium=C2 > 0.8 & EBrainMask;

    NAWMFN=[WorkingP 'ManualNAWM.nii'];
    if(exist(NAWMFN,'file'))
        CCSplenium=loadniidata(NAWMFN)>0;
    end

    RValSplenium=830;
    SparseWM=find(CCSplenium);
    RepSliIm=squeeze(MeanVol(:,:,RepSli));
    WMIm=RepSliIm+squeeze(CCSplenium(:,:,RepSli))*400;
    figure(12);clf;imagesc(mritransform(WMIm));colormap gray
    title('WM mask');
    gprint(12,[WorkingP 'WMMask.png']);
    close(12);
end
%% Partition into groups, run on every group and compare results
clear TRGroupC
TRGroupC{1}=1:nOtherFAs;
CurCGroup=1;

clear FA4DF
for i=1:nOtherFAs
    FA4DF(:,:,:,i)=loadniidata(OtherFACrgFNs{i});
end
WarningStatus=warning('off','MATLAB:rankDeficientMatrix');

nFAGroups=numel(TRGroupC);
ResNFAs=cell(nFAGroups,2);
T1Res=cell(nFAGroups,3,2);
kCoeff=ones(1,2)*6;

TRGrpName=['[' regexprep(num2str(TRGroupC{CurCGroup}),'\D*','_') ']'];
% RWorkingP=[BaseWorkingP TRGrpName filesep];
MatFN=[RWorkingP 'NFARes.mat'];
if(exist(MatFN,'file'))
    disp('Skipping');
    load(MatFN);
%     continue;
end

TRGroup=TRGroupC{CurCGroup};
Needed=TRGroup; % [1:4];
FAs=FAsF(Needed);
TRs=TRsF(Needed);
FA4D=FA4DF(:,:,:,Needed);
nFAs=numel(FAs);
SFA=size(FA4D);
ANFAs=zeros(SFA(3),nFAs);
ValMsk=all(FA4D>100,4);
% FBrainMaskS=FBrainMask & ValMsk;

%     CalcRelaxForVolNFA;
ONFAs=FAs;
%     load(MatFN,'ONFAs','FFA2DData','FBrainMaskS');
%     load(MatFN,'ONFAs');

FBrainMaskS=BrainMask;

% WMFA2DData=Reshape4d22d(FA4D,CCSplenium);

% [RelaxCG{1:3}]=CalcT1byFAfw2(WMFA2DData',FAs,TRs);
% MWMa=getKthElement(mean(RelaxCG{1}),1);
% kFA(1)=sqrt(MWMa/RValSplenium);

% [RelaxCG{1:3}]=CalcT1byFAfw2(WMFA2DData',ONFAs,TRs);
% MWMa=getKthElement(mean(RelaxCG{1}),1);
% kFA(2)=sqrt(MWMa/RValSplenium);

% FAsK=kFA(1)*FAs;
% ONFAsK=kFA(2)*ONFAs;
% save(MatFN,'ONFAs','ONFAsK','kFA','FAsK');

FFA2DData=Reshape4d22d(FA4D,FBrainMaskS>-1);

[RelaxCG{1:3}]=CalcT1byFAfw2(FFA2DData',FAs,TRs);
FRelax4D3=Reshape2DCto4D(RelaxCG,FBrainMaskS>-1);

Raw2Nii(squeeze(FRelax4D3(:,:,:,1)),[RWorkingP 'T13DOFA.nii'],'float32', MeanFN);
Raw2Nii(squeeze(FRelax4D3(:,:,:,2)),[RWorkingP 'PD3DOFA.nii'],'float32', MeanFN);
Raw2Nii(squeeze(FRelax4D3(:,:,:,3)),[RWorkingP 'RMS3DOFA.nii'],'float32', MeanFN);

if(DoDeviations)
    CalcRelaxForVolNFA;

    [RelaxCG{1:3}]=CalcT1byFAfw2(FFA2DData',ONFAs,TRs);
    FRelax4DN=Reshape2DCto4D(RelaxCG,FBrainMaskS>-1);

    Raw2Nii(squeeze(FRelax4DN(:,:,:,1)),[RWorkingP 'T13DNFA.nii'],'float32', MeanFN);
    Raw2Nii(squeeze(FRelax4DN(:,:,:,2)),[RWorkingP 'PD3DNFA.nii'],'float32', MeanFN);
    Raw2Nii(squeeze(FRelax4DN(:,:,:,3)),[RWorkingP 'RMS3DNFA.nii'],'float32', MeanFN);

    T1Res{CurCGroup,2,1}=[RWorkingP 'T13DNFA.nii'];
    T1Res{CurCGroup,2,2}=Iterative_n3_b1_clean_mask_FN(3,T1Res{CurCGroup,2,1},FMaskFN,CalcForce,[]);
end


T1Res{CurCGroup,1,1}=[RWorkingP 'T13DOFA.nii'];
T1Res{CurCGroup,1,2}=Iterative_n3_b1_clean_mask_FN(3,T1Res{CurCGroup,1,1},FMaskFN,CalcForce,[]);
kCoeff=[];
if(DoGlobal)
    T1Res{CurCGroup,1,3}=[RWorkingP 'T13DOFA_N3k.nii'];
    if(DoDeviations)
            T1Res{CurCGroup,2,3}=[RWorkingP 'T13DNFA_N3k.nii'];
    end
    for i=1:(1+DoDeviations)
        Tmp=loadniidata(T1Res{CurCGroup,i,2});
        MWMc=mean(Tmp(CCSplenium));
        kCoeff(i)=RValSplenium/MWMc;

        Raw2Nii(Tmp*kCoeff(i),T1Res{CurCGroup,i,3},'float32', MeanFN);
    end
end
save(MatFN,'ONFAs','kCoeff');

disp('DCET1_RelaxForSubjectf finished');