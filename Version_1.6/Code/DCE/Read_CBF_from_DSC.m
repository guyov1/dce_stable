% DCE map size to fit to
num_rows = 256;
num_cols = 256;

% Dimension of slices (of CBF map)
slice_dim = 3;

% Load CBF map (.nifti file)
%Map = loadniidata(...
%    'D:\users\guyn\DCE_DSC_Data\1SHAHAF_TOMRY\Study20081116_103251\SHTO_Se20_DSC-Perfusion-1 min left\Penguin Results\DSC_Results_sCBFlr.nii');

% Load MTT map (.nifti file)
Map = loadniidata(...
    'D:\users\guyn\DCE_DSC_Data\1SHAHAF_TOMRY\Study20081116_103251\SHTO_Se20_DSC-Perfusion-1 min left\Penguin Results\DSC_Results_sMTTlr.nii');

% Round the values
Map = round(Map);

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
subplot(1,2,1);
imshow(Map_norm(:,:,8));
subplot(1,2,2);
imshow(Map_resized_norm(:,:,8));

% In the example I worked on, DCE volumes 1:12 fitted DSC volumes 2:13
DSC_Map_norm_final = Map_resized_norm(:,:,2:end);

% CVI are the representing voxels indices
% For the specified example, I found CVI(12) to be a good candidate
Rep_Voxel_idx = CVI(26);

% Normalize map according to rep voxel (rep voxel MTT is '1')
DSC_Map_norm_final_rep_voxel = DSC_Map_norm_final ./ DSC_Map_norm_final(Rep_Voxel_idx);

% DEBUG
% Show a normalized map of one slice

figure;imshow(mat2gray(DSC_Map_norm_final_rep_voxel(:,:,8)))

check = squeeze(DSC_Map_norm_final_rep_voxel(:,:,8));
M = zeros([sum(sum(check>0)) 3]);
[M(:,1),M(:,2),M(:,3)] = find(check);

n = 30;
xi=linspace(min( M(:,1)),max( M(:,1)),n);
yi=linspace(min(M(:,2)),max(M(:,2)),n);
[XI, YI]=meshgrid(xi,yi);
ZI = griddata(M(:,1),M(:,2),M(:,3),XI,YI);
figure;
contourf(XI,YI,ZI);
colorbar;

