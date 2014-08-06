WorkingP='\\fmri-t9\users\Moran\OptDCEinMS\MS-IT-MTX\Sub01_ARIE_CHEN\Study20140520_102624_baseline\DCE\standard\ArCh_20140520\';
MeanFN=[WorkingP 'DCEMean.nii'];
AllData=load([WorkingP 'AfterCTC.mat']);
mixedsig=AllData.CTC2D;
WhichTimePoints=1:30;
WhichVoxels=randperm(size(mixedsig,1));
WhichVoxels=WhichVoxels(1:1000);
SmallData=mixedsig(WhichVoxels,WhichTimePoints);
nICs=10;
[Out1, Out2, Out3] = fastica(SmallData,'numOfIC', nICs);

figure(3);clf;
for i=1:nICs
    subplot(nICs,1,i);
    plot(Out1(i,:));
end
%% After looking:
ArtIdx=4; % Select component Nomber
ThreshP=0.95;
nMat=mixedsig(:,WhichTimePoints)-repmat(mean(mixedsig(:,WhichTimePoints),2),[1 numel(WhichTimePoints)]);
Norms=sqrt(sum(nMat.*nMat,2));
nMat=nMat./repmat(Norms,[1 numel(WhichTimePoints)]);
nIC=Out1(ArtIdx,:)'-mean(Out1(ArtIdx,:)');
nIC=nIC./norm(nIC);
[MxVal, MxValIdx]=max(abs(nIC));
nIC=nIC*sign(nIC(MxValIdx));
C=mixedsig(:,WhichTimePoints)*nIC;
CorrMap=AllData.Msk2*0;
CorrMap(AllData.Msk2)=C;
S=sort(C);
ThreshVal=S(floor(numel(S)*ThreshP));
BMap=CorrMap>ThreshVal;
Raw2Nii(BMap,[WorkingP 'ManualArtMask.nii'],'float32',MeanFN);