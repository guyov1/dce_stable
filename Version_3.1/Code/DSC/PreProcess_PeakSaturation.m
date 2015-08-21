function [BolusSatMat ref] = PreProcess_PeakSaturation(voxels_data_4D,Filter_struct,Mask_prev)

BolusSatMat=ones(size(Mask_prev)).*Mask_prev;
ref=[];
data_inds_mask_prev=find(Mask_prev);
voxels_data_2D=reshape(voxels_data_4D,numel(Mask_prev),size(voxels_data_4D,4));
for ii=1:length(data_inds_mask_prev)
   voxel_ind=data_inds_mask_prev(ii);
   DSC_curve_voxel=voxels_data_2D(voxel_ind,:);
   [bolus_val bolus_ind]=min(DSC_curve_voxel);
   [min_pks,min_locs]=findpeaks(-DSC_curve_voxel);
   [max_pks,max_locs]=findpeaks(DSC_curve_voxel);
   %find if there's another local min in the neiborhood of the bolus:
   neib_inds=[bolus_ind-4:bolus_ind+4];
   neib_inds(find(neib_inds==bolus_ind))=[];
   max_locs_within_min_peaks=max_locs(intersect(find(max_locs>neib_inds(1)),find(max_locs<neib_inds(end))));
   neib_min_locs=intersect(min_locs,neib_inds);
   if ~isempty(neib_min_locs) %There is another minimum --> saturation effect
       [row col slice]=ind2sub(size(BolusSatMat),voxel_ind);
%        BolusSatMat(row,col,slice)=1;
      for local_max_loc=max_locs_within_min_peaks
          if local_max_loc<bolus_ind
              local_min_loc=neib_min_locs(1);
              
          else
              local_min_loc=neib_min_locs(end);
          end
          % if the 2nd min peak is close to the 1st peak - consider as
          % saturation
          if DSC_curve_voxel(local_min_loc)-bolus_val < DSC_curve_voxel(local_max_loc)-DSC_curve_voxel(local_min_loc)
              BolusSatMat(row,col,slice)=0;
          end
      end
   end
end
