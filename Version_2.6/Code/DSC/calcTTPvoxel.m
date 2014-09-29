function [TTP_voxel_time] = calcTTPvoxel(DSC_curve,TimeBetweenSamples,BolusStartGlobal)

% Find bolus start:
Ps=zeros(1,numel(DSC_curve))+2;
% We use the t-test to get the biggest probability that the distribution of the sample
% is diffrent than the rest of the test ( -> smallest Ps value)
for i=BolusStartGlobal-3:BolusStartGlobal+3    %3:min(numel(DSC_curve)-2) %Take the minimum out of 2 minutes frame to 2 frames before the end
    [h Ps(i)]=ttest2(DSC_curve(1:i),DSC_curve((i+1):end),[],[],'unequal');
end
mLPs=-log(Ps);
% figure;plot(1:numel(MedTC),MedTC,'b',1:numel(MedTC),mLPs.*(max(MedTC)-min(MedTC))./(max(mLPs)-min(mLPs))+min(MedTC),'r')
[Tmp, BolusStart]=max(mLPs);

% Find bolus peak:
[~,BolusPeakFromStart]=min(DSC_curve(1,1,1,BolusStart:end));
BolusPeak=BolusPeakFromStart+BolusStart-1;

TTP_voxel_samples=BolusPeak-BolusStart+1;
TTP_voxel_time=TimeBetweenSamples*TTP_voxel_samples; % two secs between time samples
