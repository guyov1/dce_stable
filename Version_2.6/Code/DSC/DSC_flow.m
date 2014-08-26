% DSC maps flow

clc;
clear all;
close all;

AIF_auto=0;
%% get nii/dicom:  %%%%%%%%%%%%%%%%%

nii_path='D:\users\chenh\DSC_project\data\sample1\nii_files\sample1_nii.nii';
nii_path='C:\Users\Chen\Documents\MSc studies\TAU\MRI project\development\sample1_nii.nii';
[Out Header A]=loadniidata(nii_path);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-process 4D data - filtering and masking bad voxels:  %%%%%%
% filter to choose from
% 1. 'LowAtZero'
% 2. 'LowAverageBaseline'
% 3. 'NoBolusPeak' (a bit problematic)
% 4. 'PeakOutOfBolusRange'
% 5. 'BigFluctuations'

% Choose which filters to use from the pre-processing ('1') or not use
% ('0'):
filter_is_on=[1 0 0 0  1];
% arrange the image data in the proper orientation:
voxels_data_4D=flipdim(permute(Out,[2 1 3 4]),1);

[Mask,data4D_filtered,Filter]=data4D_mask_filtering(voxels_data_4D,filter_is_on);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% convert the filtered 4D data to concentration c(t):   %%%%%%%%%%%%
name_str=Header.hdr.hk.db_name;
TE=str2double(name_str(strfind(name_str,'TE')+3:strfind(name_str,'TE')+4));
if TE<0 || TE>1000
    error('TE has an invalid value: ',TE);
end
baseline_edges=Get_baseline_edges();
[concentration_4D] = Intens2concentration_4D(voxels_data_4D,baseline_edges,TE,Mask);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if AIF_auto
    % choose AIF automatically:
    AIF=choose_auto_AIF(concentration_4D(:,:,:,baseline_edges(1):end),Mask);
   
else
    % just for check:
    AIF=concentration_4D(83,67,8,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 N_data_voxels=length(Filter(5).data_voxels_indices);
    N_DSC_curves_to_use=200;
    check_vec=sort(random('unid',N_data_voxels,1,N_DSC_curves_to_use));
    [row_col_slice]=Filter(5).data_voxels_row_col_slice(check_vec(:),:);
    DSC_curve=zeros(1,size(voxels_data_4D,4));
    for ii=1:N_DSC_curves_to_use
        DSC_curve=DSC_curve+squeeze(voxels_data_4D(row_col_slice(ii,1),row_col_slice(ii,2),row_col_slice(ii,3),:))';
        if ii<=N_DSC_curves_to_use/(10/9)
            AIF=DSC_curve;
        end
    end
    DSC_curve=DSC_curve(baseline_edges(1):end)/N_DSC_curves_to_use;
    AIF=AIF(baseline_edges(1):end)/N_DSC_curves_to_use*(10/9);
    [Ct] =Intens2concentration_voxel(DSC_curve,baseline_edges,TE);
    AIF=Intens2concentration_voxel(AIF,baseline_edges,TE);


Rt_voxel = SVD_solve_voxel(AIF,Ct,2,'constant');
a=5;

