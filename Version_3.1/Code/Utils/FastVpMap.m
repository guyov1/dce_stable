WorkingP='\\fmri-t9\users\Moran\OptDCEinMS\MS-IT-MTX\Sub01_ARIE_CHEN\Study20140520_102624_baseline\DCE\long\ArCh_20140520\';
%%
Params=load([WorkingP 'PKM.mat']);
InterpolationFactor=ceil(Params.TimeBetweenDCEVolsFinal*Params.Options.SubSecondResolution);
HTs=0:Params.TimeBetweenDCEVolsFinal/InterpolationFactor:Params.GoodTs(end)*60;
STsI=1:InterpolationFactor:numel(HTs);
% figure;plot(HTs,Params.HAIF);hold on;plot(Params.GoodTs*60,Params.HAIF(STsI),'b*');
[Tmp, MxI]=max(Params.HAIF);
PeakIdsx=floor(MxI/InterpolationFactor+1)-1:ceil(MxI/InterpolationFactor+1)+1;
CTC4D=loadniidata([WorkingP 'CTC4D.nii']);
FastVp=max(CTC4D(:,:,:,PeakIdsx),[],4);
MeanFN=[WorkingP 'DCEMean.nii'];
Raw2Nii(FastVp,[WorkingP  'AutoArtBAT' filesep 'FastVp.nii'],'float32', MeanFN);