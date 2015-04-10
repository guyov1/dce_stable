function [voxels_peak_in_range ref] = PreProcess_PeakOutOfBolusRange(voxels_data_4D,Filter_struct,Mask_prev)

    % The output is "1" for voxels with correct location of minimums, "-1"
    % for wrong location of minimums.
    
    % Define how many minimums will be considered for decide "1" or "-1"
    num_smallest_elements_to_consider=3;
    bolus_range_from_mean=Filter_struct.bolus_range_from_mean;
    [~, mean_bolus_point_ind] = min(Filter_struct.mean_DSC_curve);
    range_bolus_search = [mean_bolus_point_ind-bolus_range_from_mean : mean_bolus_point_ind+bolus_range_from_mean];

    [N_rows,N_cols,N_slices,N_time_points]=size(voxels_data_4D);
    voxels_peak_in_range=ones(N_rows,N_cols,N_slices);
    
    for row=1:N_rows
        for col=1:N_cols
            for slice=1:N_slices
                voxel_time_vec_sorted=sort(voxels_data_4D(row,col,slice,:),4); %sort along the time curve
                % find the smallest elements along the time sample:
                smallest_elements_inds=find(voxels_data_4D(row,col,slice,:)<=voxel_time_vec_sorted(num_smallest_elements_to_consider));
                % check if any of the smallest elements is outside a
                % defined "bolus range"
                if length(intersect(smallest_elements_inds,range_bolus_search))<length(smallest_elements_inds)
                    voxels_peak_in_range(row,col,slice)=-1;
                end
            end
        end
    end

    
    ref=[];
end