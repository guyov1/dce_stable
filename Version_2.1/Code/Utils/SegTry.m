% BaseP='C:\STRDCE\John\Database\DCEOut\BeDa_20130617\';
BaseP='C:\STRDCE\John\Database\DCEOut\HaMo_20130714\';
ASegP=[BaseP 'Seg_BaselineNoBadSli2' filesep];
MeanFN=[BaseP 'DCEMean.nii'];
FMaskFN=[BaseP 'FBrainMsk.nii'];
A=loadniidata([ASegP 'ForSeg.nii']);
M=loadniidata(FMaskFN)>0;
M=A>0;
S=A;
GoodSlicesS=find(isfinite(squeeze(gsum(M,1:2))'));
se=strel('disk',4,8);
H = fspecial('disk',10);
for i=GoodSlicesS
    S(:,:,i)=imfilter(squeeze(A(:,:,i)),H,'replicate');
end
N=median(A(A>0))*(A./S);
N(~M)=NaN;

TempFN=[BaseP 'Temp.nii'];
Raw2Nii(N,TempFN,'float32', MeanFN);
DCEMeanSegP2=SPM_Segment(TempFN,true,[],false);

mricronx({TempFN,[BaseP 'Seg_Temp\c1ForSeg.nii'],[BaseP 'Seg_Temp\c2ForSeg.nii'],[BaseP 'Seg_Temp\c3ForSeg.nii']},[0 1 2 3]);