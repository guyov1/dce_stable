function [max_min_in_bolus_range_diff ref] = PreProcess_WrongBolusPeak(voxels_data_4D,Filter_struct,Mask_prev)
    
    ref=[];

    first_baseline_point_for_mask=Filter_struct.first_baseline_point_to_include;
    last_baseline_point_for_mask=Filter_struct.first_baseline_point_to_include+Filter_struct.baseline_length_for_mask-1;
    last_bolus_range_point_for_mask=last_baseline_point_for_mask+Filter_struct.bolus_range_length-1;
%     [voxels_min_vals_in_time voxels_min_indices_in_time] = min(voxels_data_4D,[],4);
    voxels_baseline_mean = mean(voxels_data_4D(:,:,:,first_baseline_point_for_mask:last_baseline_point_for_mask),4);
    voxels_bl_diff=voxels_data_4D-repmat(voxels_baseline_mean,[1 1 1 size(voxels_data_4D,4)]);
    [max_vals max_inds]=max(voxels_data_4D(:,:,:,last_baseline_point_for_mask+1:last_bolus_range_point_for_mask),[],4);
    [min_vals min_inds]=min(voxels_data_4D(:,:,:,last_baseline_point_for_mask+1:last_bolus_range_point_for_mask),[],4);    
    min_vals_concent=-log(max_vals./voxels_baseline_mean); %not normalized - doesn't matter for filtering. 
    max_vals_concent=-log(min_vals./voxels_baseline_mean); %not normalized - doesn't matter for filtering    
    max_min_in_bolus_range_diff=max_vals_concent-abs(min_vals_concent);
%     voxels_abs=abs(voxels_bl_diff);
%     [max_abs_vals inds]=max(voxels_abs,[],4);
%     inds_absolute=inds+last_baseline_point_for_mask;
    
% 
%     % ratio between min point (peak of bolus) and baseline:
%     % add "epsilon" to the denominator to avoid dividing 0:0 in the out-of-mask
%     % voxels:f
%     voxels_bl_min_ratio=abs(voxels_baseline_mean./(voxels_min_vals_in_time+eps));

%     max_value_in_bolus_range=5;
end