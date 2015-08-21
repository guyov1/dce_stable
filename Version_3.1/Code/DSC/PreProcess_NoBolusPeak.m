function [voxels_bl_min_ratio ref] = PreProcess_NoBolusPeak(voxels_data_4D,Filter_struct,Mask_prev)
    
    bolus_range_from_mean=Filter_struct.bolus_range_from_mean; % num of time points to look for bolus peak from each side of the bolus peak mean over all voxels
    first_baseline_point_for_mask=Filter_struct.first_baseline_point_to_include;
    last_baseline_point_for_mask=Filter_struct.first_baseline_point_to_include+Filter_struct.baseline_length_for_mask-1;
    
    [~, mean_bolus_point_ind] = min(Filter_struct.mean_DSC_curve);
    range_bolus_search = [mean_bolus_point_ind-bolus_range_from_mean : mean_bolus_point_ind+bolus_range_from_mean];

    [voxels_min_vals_in_time voxels_min_indices_in_time] = min(voxels_data_4D(:,:,:,range_bolus_search),[],4);
    voxels_baseline_mean = mean(voxels_data_4D(:,:,:,first_baseline_point_for_mask:last_baseline_point_for_mask),4);

    % ratio between min point (peak of bolus) and baseline:
    % add "epsilon" to the denominator to avoid dividing 0:0 in the out-of-mask
    % voxels:
    voxels_bl_min_ratio=abs(voxels_baseline_mean./(voxels_min_vals_in_time+eps));

    
    ref=[];
end