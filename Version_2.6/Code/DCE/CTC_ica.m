WorkingP='D:\DCE\DataBase\C001_FoNa_20100603\DCE\FoNa_20100603\';
WhichTimePoints=1:30;   % What time points to use for the ICA (just to use a bit less memory)
nICs=10;                % How many Independent Components to compute (arteries should be one of the first components)
NoOfVoxelsToTakeForICA=1000; % How many random voxels to sample for ICA (just to use a bit less memory)
SlicesToIgnoreFromTop=1;    % Don't take voxels from these slices at the top
SlicesToIgnoreFromBottom=1; % Don't take voxels from these slices at the bottom

MeanFN=[WorkingP 'DCEMean.nii'];
AllData=load([WorkingP 'AfterCTC.mat']);
mixedsig=AllData.CTC2D;

IdxMap=AllData.Msk2*0;
IdxMap(AllData.Msk2)=1:size(AllData.CTC2D,1);
SlicesWithVoxels=find(squeeze(sum(sum(IdxMap,1),2))>10);
GoodSlices=(min(SlicesWithVoxels)+SlicesToIgnoreFromBottom):(max(SlicesWithVoxels)-SlicesToIgnoreFromTop);

IdxInGoodSlices=IdxMap(:,:,GoodSlices);
IdxInGoodSlices=IdxInGoodSlices(IdxInGoodSlices>0);

WhichVoxels=getKrandomSamples(IdxInGoodSlices,NoOfVoxelsToTakeForICA);
SmallData=mixedsig(WhichVoxels,WhichTimePoints);
[Out1, Out2, Out3] = fastica(SmallData,'numOfIC', nICs);

figure(3);clf;
for i=1:nICs
    subplot(nICs,1,i);
    plot(Out1(i,:));
end
%% After looking:
ArtIdx=2;       % Select component Number for arteries
ThreshP=0.95;   % How correlative should it be to that component to be included
nMat=mixedsig(:,WhichTimePoints)-repmat(mean(mixedsig(:,WhichTimePoints),2),[1 numel(WhichTimePoints)]);
Norms=sqrt(sum(nMat.*nMat,2));
nMat=nMat./repmat(Norms,[1 numel(WhichTimePoints)]);
nIC=Out1(ArtIdx,:)'-mean(Out1(ArtIdx,:)');
nIC=nIC./norm(nIC);
[MxVal, MxValIdx]=max(abs(nIC));
nIC=nIC*sign(nIC(MxValIdx));
C=nMat*nIC;
CorrMap=AllData.Msk2*0;
Nifits(AllData.Msk2)=C;
S=sort(C);
ThreshVal=S(floor(numel(S)*ThreshP));
BMap=CorrMap>ThreshVal;
Raw2Nii(BMap,[WorkingP 'ManualArtMask.nii'],'float32',MeanFN);
%% Output Nifits
Corr4D=zeros([size(CorrMap,1) size(CorrMap,2) size(CorrMap,3) nICs]);
for i=1:nICs
    nIC=Out1(i,:)'-mean(Out1(i,:)');
    nIC=nIC./norm(nIC);
    [MxVal, MxValIdx]=max(abs(nIC));
    nIC=nIC*sign(nIC(MxValIdx));
    C=nMat*nIC;
    CorrMap(AllData.Msk2)=C;
    Corr4D(:,:,:,i)=CorrMap;
end
Raw2Nii(Corr4D,[WorkingP 'IComponents.nii'],'float32',MeanFN); 