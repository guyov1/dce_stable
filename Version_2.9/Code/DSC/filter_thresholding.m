function [data_voxels_row_col_slice,data_voxels_indices,noise_voxels_indices,Mask] = filter_thresholding(data_for_filter,th,Mask_prev,ref)

    % "ref" is an extra parameter that some of the filters need.
    % in the case "ref" exist, the threshold is only a factor to multiply
    % the reference value "ref"
    if ~isempty(ref)
       th=th*ref; 
    end
    
    noise_voxels_indices = find(data_for_filter<th);
    Mask=Mask_prev;
    Mask(noise_voxels_indices)=0;
    
    [row_data,col_data,slice_data]=ind2sub(size(data_for_filter),find(Mask==1));
    data_voxels_row_col_slice=[row_data(:) col_data(:) slice_data(:)];
    data_voxels_indices = find(Mask==1);

    
%     [row_data,col_data,slice_data]=ind2sub(size(data_for_filter),find(data_for_filter>=th));
%     data_voxels_row_col_slice=[row_data(:) col_data(:) slice_data(:)];
%     data_voxels_indices=find(data_for_filter>=th);
%     
%     
%     
%     noise_voxels_indices= data_for_filter<th;
%     
%     Mask=Mask_prev;
%     Mask(noise_voxels_indices)=0;

end