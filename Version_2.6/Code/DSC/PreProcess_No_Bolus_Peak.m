function [voxels_bl_min_ratio] = PreProcess_NoBolusPeak(voxels_data_4D,Filter_struct)
    
    first_baseline_point_for_mask=Filter_struct.first_baseline_point_to_include;
    last_baseline_point_for_mask=Filter_struct.first_baseline_point_to_include+Filter_struct.baseline_length_for_mask-1;

    [voxels_min_vals_in_time voxels_min_indices_in_time] = min(voxels_data_4D,[],4);
    voxels_baseline_mean = mean(voxels_data_4D(:,:,:,first_baseline_point_for_mask:last_baseline_point_for_mask),4);

    % ratio between min point (peak of bolus) and baseline:
    % add "epsilon" to the denominator to avoid dividing 0:0 in the out-of-mask
    % voxels:
    voxels_bl_min_ratio=abs(voxels_baseline_mean./(voxels_min_vals_in_time+eps));

end