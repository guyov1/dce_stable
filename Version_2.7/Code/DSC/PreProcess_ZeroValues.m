function [zero_lengths ref]=PreProcess_ZeroValues(voxels_data_4D,Filter_data,Mask_prev)

% We don't want voxels whose intensity curves contain a few
% consecutive zeros:
voxels_data_4D_permute=permute(voxels_data_4D,[4,3,2,1]);
[row_data col_data slice_data time_data]=ind2sub(size(voxels_data_4D_permute),find(voxels_data_4D_permute==0));
full_data=[row_data col_data slice_data time_data];


zero_lengths=1;
ref=1;