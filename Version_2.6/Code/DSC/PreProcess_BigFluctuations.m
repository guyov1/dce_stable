function [voxels_mean_minus_abs_diff_in_time , total_mean_minus_abs_diff_in_time] = PreProcess_BigFluctuations(voxels_data_4D,Filter_struct,Mask_prev)

    [N_rows,N_cols,N_slices,N_time_points]=size(voxels_data_4D);
    
    voxels_3D_diff_in_time=zeros(N_rows,N_cols,N_slices,N_time_points-1);
    for t=1:N_time_points-1
        voxels_3D_diff_in_time(:,:,:,t)=voxels_data_4D(:,:,:,t+1)-voxels_data_4D(:,:,:,t);
    end
    voxels_3D_abs_diff_in_time=abs(voxels_3D_diff_in_time);
    voxels_vec_abs_diff_in_time=voxels_3D_abs_diff_in_time(:);
    voxels_mean_minus_abs_diff_in_time=-mean(voxels_3D_abs_diff_in_time,4);
    
    total_mean_minus_abs_diff_in_time=-mean(voxels_vec_abs_diff_in_time(Mask_prev>0));
end