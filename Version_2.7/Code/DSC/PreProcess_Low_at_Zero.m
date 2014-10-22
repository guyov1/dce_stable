function [voxels_first_time_sample] = PreProcess_LowAtZero(voxels_data_4D,Filter_data)

% Thresholding to delete blank voxels:
% We separate voxels to "data" (over th) and "noise" (below th)
% Based on the first time sample.
    voxels_first_time_sample=voxels_data_4D(:,:,:,1);
end