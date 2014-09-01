% WorkingP='\\fmri-t9\users\Moran\OptDCEinMS\MS-IT-MTX\Sub01_ARIE_CHEN\Study20140520_102624_baseline\DCE\long\ArCh_20140520\';
WorkingP='\\fmri-t9\users\Moran\OptDCE\Healthy\03_Ifat_Harris\DCE\IfHa_20101230\';
DataP=[WorkingP 'AutoArtBAT' filesep];
PercentToUse=0.9;
%%
CTC4D=loadniidata([WorkingP 'CTC4D.nii']);
VeinsROI=loadniidata([WorkingP 'veins.nii']);
CTC2D=Reshape4d22d(CTC4D,VeinsROI>0);
load([WorkingP 'PKM.mat']);
load([WorkingP 'Params.mat']);

SampleTs=GoodTs;
TimeVec=SampleTs;
% As before, set all relevant paramters
TimeBetweenDCEVolsMin=TimeBetweenDCEVolsFinal/60;
InterpolationFactor=ceil(TimeBetweenDCEVolsFinal);

HInterpolationFactor = ceil(InterpolationFactor*Options.SubSecondResolution);
% HInterpolationFactor=ceil(InterpolationFactor*Options.SubSecRes);
Hdt                  = TimeBetweenDCEVolsMin/HInterpolationFactor;
HSampleTs            = 0:Hdt:SampleTs(end);
DTimeVec=diff(TimeVec);
F=find(DTimeVec>DTimeVec(1)*2,1);
nNormalTs=numel(TimeVec);
if(~isempty(F))
    nNormalTs=F(1);
end
HSampleTs=[0:Hdt:SampleTs(nNormalTs) TimeVec((nNormalTs+1):end)];

InspectedAIFParamsTimeFN=[WorkingP 'InspectedAIFParamsTime.mat'];
if(exist(InspectedAIFParamsTimeFN,'file')) 
    Tmp1=load(InspectedAIFParamsTimeFN);
    OldTimeVec=Tmp1.InspectedParamsTimeSamples;
    AIF_Parker9t=@(x,t) AIF_Parkerg3( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8),x(9),max(OldTimeVec))*x(2);
else
    AIF_Parker9t=@(x,t) AIF_Parkerg3( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8),x(9),max(TimeVec))*x(2);
end

% Create Parker's AIF
HAIF=AIF_Parker9t(OutAIFParam,HSampleTs);

DF=diff(GoodTs);
nMainDCE=find(DF>DF(1)*1.1,1);
nTimePointsToUse=numel(SampleTs);

IntAIF=trapz(HSampleTs,HAIF);
IntVeins=trapz(SampleTs,CTC2D(:,1:nTimePointsToUse)');
NVals=numel(IntVeins);
Idx=floor(PercentToUse*NVals);
Sorted=sort(IntVeins);
MxIntVeins=Sorted(Idx);
NewFactor=MxIntVeins/IntAIF;
NewHAIF=HAIF*NewFactor;

[Mx, MxI]=max(HAIF);
HIs=(MxI-HInterpolationFactor):(MxI+HInterpolationFactor);
IntAIFBol=trapz(HSampleTs(HIs),HAIF(HIs));
for i=1:size(CTC2D,1)
    CurCTC=CTC2D(i,:);
    [Mx, MxI]=max(CurCTC);
    Is=(MxI-1):(MxI+1);
    IntVeinsBol(i)=trapz(SampleTs(Is),CTC2D(i,Is)');
end
NVals=numel(IntVeinsBol);
Idx=floor(PercentToUse*NVals);
Sorted=sort(IntVeinsBol);
MxIntVeinsBol=Sorted(Idx);
NewFactorBol=MxIntVeinsBol/IntAIFBol;
NewHAIFBol=HAIF*NewFactorBol;

FromJimCoeff=NewFactor/AIFAmpCoeff;
FromJimCoeffBol=NewFactorBol/AIFAmpCoeff;

figure;
plot(SampleTs,CTC2D(:,1:nTimePointsToUse)');
hold on;
h(1)=plot(HSampleTs,HAIF*AIFAmpCoeff,'b','LineWidth',3);
h(2)=plot(HSampleTs,NewHAIF,'m','LineWidth',3);
h(3)=plot(HSampleTs,NewHAIFBol,'r','LineWidth',3);
title(['New factor: ' num2str(FromJimCoeff) ', ' num2str(FromJimCoeffBol)]);
legend(h,{'Jim','Sorbronne','Sorbronne Bol'});


MeanFN=[WorkingP 'DCEMean.nii'];

TmpA=loadniidata([DataP 'KtransFinalN.nii']);
Raw2Nii(TmpA/FromJimCoeffBol,[DataP 'KtransFinalNS.nii'],'float32', MeanFN);
TmpA=loadniidata([DataP 'VpFinalN.nii']);
Raw2Nii(TmpA/FromJimCoeffBol,[DataP 'VpFinalNS.nii'],'float32', MeanFN);

TmpA=loadniidata([DataP 'KtransFinalN.nii']);
Raw2Nii(TmpA/FromJimCoeff,[DataP 'KtransFinalNSB.nii'],'float32', MeanFN);
TmpA=loadniidata([DataP 'VpFinalN.nii']);
Raw2Nii(TmpA/FromJimCoeff,[DataP 'VpFinalNSB.nii'],'float32', MeanFN);

% Save AIF factor
save([WorkingP 'AIF_Sourbron_Factor.mat'],'NewFactor')