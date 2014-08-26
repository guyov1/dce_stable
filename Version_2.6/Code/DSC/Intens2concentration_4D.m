function [concentration_4D Mask] = Intens2concentration_4D(voxels_data4D,baseline_edges,last_sample,TE,Mask,handles)
    
    [N_rows,N_cols,N_slices,N_time_points]=size(voxels_data4D);
    data_voxels_row_col_slice=handles.data_voxels_row_col_slice;
    concentration_4D=zeros(size(voxels_data4D));
    
    for ii=1:size(data_voxels_row_col_slice,1)
       row=data_voxels_row_col_slice(ii,1);
       col=data_voxels_row_col_slice(ii,2);
       slice=data_voxels_row_col_slice(ii,3);
       intensity_curve=voxels_data4D(row,col,slice,:);
       [Ct_curve Keep_in_Mask]=Intens2concentration_voxel(intensity_curve,baseline_edges,TE);
       if Keep_in_Mask
           concentration_4D(row,col,slice,:)=Ct_curve;
       else
           Mask(row,col,slice)=0;
           data_voxels_row_col_slice(ii,:)=[];
       end
    end
    concentration_4D=concentration_4D(:,:,:,baseline_edges(1):last_sample);
    handles.data_voxels_row_col_slice=data_voxels_row_col_slice;
%     for row=1:N_rows
%         for col=1:N_cols
%             for slice=1:N_slices
%                 if Mask(row,col,slice)>0
%                     intensity_curve=voxels_data4D(row,col,slice,:);
%                     [Ct_curve Keep_in_Mask]=Intens2concentration_voxel(intensity_curve,baseline_edges,TE);
%                     if Keep_in_Mask
%                         concentration_4D(row,col,slice,:)=Ct_curve;
%                     else
%                         Mask(row,col,slice)=0;
%                     end
%                 else
%                     concentration_4D(row,col,slice,:)=0;
%                 end
%             end
%         end
%     end


end