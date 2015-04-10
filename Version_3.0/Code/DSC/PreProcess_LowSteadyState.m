function [steady_state_baseline_diff ref] = PreProcess_LowSteadyState(voxels_data_4D,Filter_struct,Mask_prev)

ref=[];

first_baseline_point_for_mask=Filter_struct.first_baseline_point_to_include;
last_baseline_point_for_mask=Filter_struct.first_baseline_point_to_include+Filter_struct.baseline_length_for_mask-1;
voxels_baseline_mean = mean(voxels_data_4D(:,:,:,first_baseline_point_for_mask:last_baseline_point_for_mask),4);

last_steady_state_point_for_mask=Filter_struct.last_point_to_include;
first_steady_state_point_for_mask=Filter_struct.last_point_to_include-Filter_struct.steady_state_part_length+1;
voxels_steady_state_mean = mean(voxels_data_4D(:,:,:,first_steady_state_point_for_mask:last_steady_state_point_for_mask),4);


steady_state_baseline_diff = voxels_baseline_mean - voxels_steady_state_mean;