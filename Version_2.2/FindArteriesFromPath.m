if(UnderSampling==1)
    USStr='';
else
    USStr=['_' num2str(UnderSampling)];
end
WorkingPx=WorkingP;
CTCFN=[WorkingP 'AfterCTC' USStr '.mat'];
disp('Loading..');
load(CTCFN);
WorkingP=WorkingPx;
PrepareFN=[WorkingP 'AfterPrepare4D.mat'];
MeanFN=[WorkingP 'DCEMean.nii'];
if(~isfield(Options,'TimeDelayToMaskVeins'))
    Options.TimeDelayToMaskVeins=-0.5;
end
% Options=handles.Options;

disp('DCET1_PKf..');
DCET1_PKf;

disp('AIFFinderTry..');
ROIStr='Full';
DebugSuffix='_ForDB';
ArtFinder;