close all;

%% General parameters
% DCE map size to fit to
num_rows = 256;
num_cols = 256;

% Dimension of slices (of CBF/MTT map)
slice_dim = 3;

% Map to choose (CBF/MTT etc.)
Map_To_Choose = 'MTT';

% Maximum value from which we zero voxel's value
Max_allowed_val = 15;

% Slice to display in debugging
Slice_num_to_disp = 8;

%%

% Read the needed map
switch Map_To_Choose
    case 'MTT'
        % Load MTT map (.nifti file)
        Map = loadniidata(...
            'D:\users\guyn\DCE_DSC_Data\1SHAHAF_TOMRY\Study20081116_103251\SHTO_Se20_DSC-Perfusion-1 min left\Penguin Results\DSC_Results_sMTTlr.nii');
    case 'CBF'
        % Load CBF map (.nifti file)
        Map = loadniidata(...
            'D:\users\guyn\DCE_DSC_Data\1SHAHAF_TOMRY\Study20081116_103251\SHTO_Se20_DSC-Perfusion-1 min left\Penguin Results\DSC_Results_sCBFlr.nii');
    otherwise
        error('-E- Cant recognize map to choose!')
end

% Round the values
%Map = round(Map);

% Zero all voxels with value greater than 15 seconds
% I assume those voxels don't go back to baseline -> do not represent AIF
display(sprintf('-I- Removing voxels with value greater than %d',Max_allowed_val));
Map( Map(:,:,:) > Max_allowed_val ) = 0;

% Interpolate CBF map to needed size
orig_num_slices = size(Map,slice_dim);
orig_num_rows   = size(Map,1);
orig_num_cols   = size(Map,2);

% Interpolation ratio between original to target resolution
interp_ratio_rows = num_rows/orig_num_rows;
interp_ratio_cols = num_cols/orig_num_cols;

[xi,yi,zi] = meshgrid(1:num_rows,1:num_cols,1:orig_num_slices);

% New number of rows and columns
%interp2(1:num_rows,1:num_cols,CBF_Map(:,:,4),1:interp_ratio_rows:num_rows,1:interp_ratio_cols:num_cols);
%CBF_Map_resized = imresize(mat2gray(CBF_Map), [num_rows num_cols orig_num_slices])
%CBF_Map_resized = interp2(1:interp_ratio_rows:num_rows,1:interp_ratio_cols:num_cols,CBF_Map(:,:,4),(1:num_rows)',(1:num_cols),'nearest');
Map_resized = imresize(mat2gray(Map(:,:,:)), [num_rows num_cols],'nearest');

% Change to normalized (0-1) values
Map_norm         = mat2gray(Map);
Map_resized_norm = mat2gray(Map_resized);

% Display original and resized MTT/CBF map
figure;
subplot(2,2,1);
imshow(Map_norm(:,:,Slice_num_to_disp));
title(sprintf('%s normalized map, before interpolation',Map_To_Choose));
subplot(2,2,2);
imshow(Map_resized_norm(:,:,Slice_num_to_disp));
title(sprintf('%s normalized map, after interpolation',Map_To_Choose));

% In the example I worked on, DCE volumes 1:12 fitted DSC volumes 2:13
DSC_Map_norm_final = Map_resized_norm(:,:,2:end);

% CVI are the representing voxels indices
% For the specified example, I found CVI(1) to be a good candidate
%Rep_Voxel_idx = CVI(14);


%  Normalize map according to median value (ignoring 0 values)
Median_value = median(median(median(DSC_Map_norm_final(DSC_Map_norm_final>0))));
DSC_Map_norm_final_rep_voxel = DSC_Map_norm_final ./ Median_value;

% Rep voxel index will be a voxel with a median value
diff_from_median = min(abs(DSC_Map_norm_final(DSC_Map_norm_final-Median_value > 0) -Median_value));
Rep_Voxel_idx = find(DSC_Map_norm_final == (Median_value + diff_from_median) );

% Normalize map according to rep voxel (rep voxel MTT is '1')
%DSC_Map_norm_final_rep_voxel = DSC_Map_norm_final ./ DSC_Map_norm_final(Rep_Voxel_idx);

slice_idx = ceil(Rep_Voxel_idx / ( size(DSC_Map_norm_final,1) * size(DSC_Map_norm_final,2) ));
col_idx   = ceil( ( Rep_Voxel_idx - (slice_idx-1)*(65536) ) / size(DSC_Map_norm_final,1));
row_idx   = mod( int64(Rep_Voxel_idx - (slice_idx-1)*(65536)) ,size(DSC_Map_norm_final,1));

% DEBUG
% Show a normalized map of one slice

%figure;
subplot(2,2,3);

imshow(mat2gray(DSC_Map_norm_final_rep_voxel(:,:,Slice_num_to_disp)));
hold on;
plot(row_idx,col_idx,'+');
%title(sprintf('%s map, after normalizing to rep voxel with norm value: %f',Map_To_Choose,DSC_Map_norm_final(Rep_Voxel_idx)));
title(sprintf('%s map, after normalizing to median value: %f',Map_To_Choose,Median_value));
hold off;

check = squeeze(DSC_Map_norm_final_rep_voxel(:,:,Slice_num_to_disp));
M = zeros([sum(sum(check>0)) 3]);
[M(:,1),M(:,2),M(:,3)] = find(check);

n = 30;
xi=linspace(min( M(:,1)),max( M(:,1)),n);
yi=linspace(min(M(:,2)),max(M(:,2)),n);
[XI, YI]=meshgrid(xi,yi);
ZI = griddata(M(:,1),M(:,2),M(:,3),XI,YI);
%figure;
subplot(2,2,4);
contourf(XI,YI,ZI);
colorbar;

