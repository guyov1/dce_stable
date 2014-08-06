WorkingP='\\fmri-t9\users\Moran\OptDCEinMS\MS-IT-MTX\Sub01_ARIE_CHEN\Study20140520_102624_baseline\DCE\long\ArCh_20140520\';
%%
CTC4D=loadniidata([WorkingP 'CTC4D.nii']);
VeinsROI=loadniidata([WorkingP 'Veins.nii']);
CTC2D=Reshape4d22d(CTC4D,VeinsROI>0);
Params=load([WorkingP 'PKM.mat']);
GeneralDataFN=[WorkingP 'Params.mat'];
Params2=load(GeneralDataFN);
InterpolationFactor=ceil(Params.TimeBetweenDCEVolsFinal*Params.Options.SubSecondResolution);
HTs=0:Params.TimeBetweenDCEVolsFinal/InterpolationFactor:Params.GoodTs(end)*60;
STsI=1:InterpolationFactor:numel(HTs);
HAIF=Params.HAIF*Params2.AIFAmpCoeff;
IntAIF=trapz(HAIF)/InterpolationFactor;
MxIntVeins=max(trapz(CTC2D'));
NewFactor=MxIntVeins/IntAIF;
NewHAIF=HAIF*NewFactor;
figure;plot(Params.GoodTs*60,CTC2D');hold on;plot(HTs,HAIF,'b','LineWidth',3);plot(Params.GoodTs*60,HAIF(STsI),'b*');
plot(HTs,NewHAIF,'m','LineWidth',3);
title(['New factor: ' num2str(NewFactor)]);