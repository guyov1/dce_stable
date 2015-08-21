function [TTP] = calcTTP(handles,TimeBetweenSamples)

TTP=zeros(size(handles.Mask));
data4D_init=handles.data4D_init;
data_voxels_row_col_slice=handles.data_voxels_row_col_slice;

for ii=1:size(data_voxels_row_col_slice,1)
    row=data_voxels_row_col_slice(ii,1);
    col=data_voxels_row_col_slice(ii,2);
    slice=data_voxels_row_col_slice(ii,3);
    DSC_curve=data4D_init(row,col,slice,:);
    
    % Find bolus start:
    Ps=zeros(1,numel(DSC_curve))+2;
    % We use the t-test to get the biggest probability that the distribution of the sample
    % is diffrent than the rest of the test ( -> smallest Ps value)
    for i=3:min(numel(DSC_curve)-2) %Take the minimum out of 2 minutes frame to 2 frames before the end
        [h Ps(i)]=ttest2(DSC_curve(1:i),DSC_curve((i+1):end),[],[],'unequal');
    end
    mLPs=-log(Ps);
    % figure;plot(1:numel(MedTC),MedTC,'b',1:numel(MedTC),mLPs.*(max(MedTC)-min(MedTC))./(max(mLPs)-min(mLPs))+min(MedTC),'r')
    [Tmp, BolusStart]=max(mLPs);
    
    % Find bolus peak:
    BolusPeak=min(DSC_curve);
    
    TTP_voxel_samples=BolusPeak-BolusStart+1;
    TTP_voxel_time=TimeBetweenSamples*TTP_voxel_samples; % two secs between time samples
    TTP(row,col,slice)=TTP_voxel_time;
end