function [voxels_baseline_mean ref] = PreProcess_LowAverageBaseline(voxels_data_4D,Filter_struct,Mask_prev)
        
    first_baseline_point_for_mask=Filter_struct.first_baseline_point_to_include;
    last_baseline_point_for_mask=Filter_struct.first_baseline_point_to_include+Filter_struct.baseline_length_for_mask-1;
    
    voxels_baseline_mean = mean(voxels_data_4D(:,:,:,first_baseline_point_for_mask:last_baseline_point_for_mask),4);

    ref=[];
end