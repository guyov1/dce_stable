%% read nii and create curves for each voxel.
clear all;
close all;
nii_path='D:\users\chenh\DSC_project\data\sample1\nii_files\sample1_nii.nii';
[Out Header A]=loadniidata(nii_path);

% arrange the image data in the proper orientation:
voxels_data_4D=flipdim(permute(Out,[2 1 3 4]),1);
[N_rows,N_cols,N_slices,N_time_points]=size(voxels_data_4D);

first_phase_included=10;
last_phase_included=92;


%% Attempt to plot DSC curve for example
col=47;
row=48;
slice=7;

Ct=voxels_data_4D(col,row,slice,first_phase_included:last_phase_included);
Ct=Ct(:);
% figure;plot(Ct);
%% check average of Ct over all voxels in a slice
slice=4;
% build a mask for voxels which has non-zero DSC curve:
non_zero_Ct_voxels=ones(N_rows,N_cols,N_slices);
for ii=1:N_rows
    for jj=1:N_cols
        if voxels_data_4D(ii,jj,slice,:)== zeros(1,1,1,N_time_points)
            non_zero_Ct_voxels(ii,jj,slice)=0;
        end
    end
end
% self check:
% figure;
% for ii=1:N_slices
%    subplot(3,5,ii);
%    imshow(non_zero_Ct_voxels(:,:,ii));
% end

%%%%%%%%  FIGURES of some DSC curves
% figure;
% for ii=1:20:121
%     for jj=1:20:121
%         DSC_curve=voxels_data_4D(ii,jj,13,:);
%         DSC_curve=DSC_curve(:);
%         plot(DSC_curve);
%         hold all;
%     end
% end
% hold off;
%%%%%%%%%%%%%%%%%
%% statistics
% histogram of intensities of first time samples, for all voxels:
figure;
% hist(voxels_data_4D(:,:,:,1:5),100);
voxels_data_vec=voxels_data_4D(:);
voxels_data_vec_first=reshape(voxels_data_4D(:,:,:,1:1),1,numel(voxels_data_4D(:,:,:,1:1)));
voxels_4D_mean=mean(voxels_data_vec_first);
voxels_4D_stdev=std(voxels_data_vec_first);

%% MASK1 - thresholding to delete blank voxels:
% We separate voxels to "data1" (over th) and "noise" (below th)
% Mask1 is based on the first time sample.
th1=500; % th for absolute value of signal intensity
[row_data1,col_data1,slice_data1]=ind2sub(size(voxels_data_4D(:,:,:,1)),find(voxels_data_4D(:,:,:,1)>=th1));
[row_noise,col_noise,slice_noise]=ind2sub(size(voxels_data_4D(:,:,:,1)),find(voxels_data_4D(:,:,:,1)<th1));
data1_indices=find(voxels_data_4D(:,:,:,1)>=th1);
noise_indices=find(voxels_data_4D(:,:,:,1)<th1);
N_data1_voxels=length(data1_indices)

mask1=ones(N_rows,N_cols,N_slices);
mask1(noise_indices)=0;

%%%%%%%% MASK1 FIGURES:
% figure;
% for ii=1:N_slices
%     subplot(3,5,ii);
%     imshow(mask1(:,:,ii));
%     title(['slice=',num2str(ii)]);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%

% applying Mask1 for the full 4D data:
voxels_data_4D_mask1=zeros(N_rows,N_cols,N_slices,N_time_points);
for t=1:N_time_points
   voxels_data_4D_mask1(:,:,:,t)=voxels_data_4D(:,:,:,t).*mask1;
end
% voxels_data_4D_mask1(noise_indices)=0;

% statistics for the data voxels after applying Mask1:
mean_DSC_curve=zeros(1,N_time_points);
std_DSC_curve=zeros(1,N_time_points);
for t=1:N_time_points
    voxels_data_3D_mask1_t=voxels_data_4D_mask1(:,:,:,t);
    voxels_data_3D_mask1_t_vec=voxels_data_3D_mask1_t(find(voxels_data_3D_mask1_t>0));
    mean_DSC_curve(t)=mean(voxels_data_3D_mask1_t_vec);
    std_DSC_curve(t)=std(voxels_data_3D_mask1_t_vec);
%     mean_DSC_curve(t)=sum(sum(sum(voxels_data_4D_mask1(:,:,:,t))))/N_data1_voxels;
%     mean_DSC_curve(t)=mean(voxels_data_4D(data1_indices));
end
%plot some random DSC curves and the mean of ALL DSC curves
N_DSC_curves_to_show=30;
check_vec=sort(random('unid',N_data1_voxels,1,N_DSC_curves_to_show));

if 0
    figure;
    for ii=1:length(check_vec)
        check_ind=check_vec(ii);
        plot(squeeze(voxels_data_4D_mask1(row_data1(check_ind),col_data1(check_ind),slice_data1(check_ind),:)));
        hold all;
    end
    std_part_to_show=1;
    title([num2str(N_DSC_curves_to_show),' random DSC curves, of voxels passed mask1. blue = mean of passed voxels + ',num2str(std_part_to_show),'std. green = mean of ALL voxel before Mask1']);
    errorbar(mean_DSC_curve,std_part_to_show*std_DSC_curve,'--b','LineWidth',2);
    % plot(mean_DSC_curve,'--b','LineWidth',3);
    plot(mean_DSC_curve_total,'-.g','LineWidth',3);
    hold off;
end


%% MASK 2 - take out voxels without bolus arrival:
first_baseline_point_to_include=first_phase_included;
baseline_length_for_mask2=10;
last_baseline_point_mask2=first_baseline_point_to_include+baseline_length_for_mask2;

% ratio between min point (peak of bolus) and baseline:
[voxels_min_vals_in_time voxels_min_indices_in_time] = min(voxels_data_4D_mask1,[],4);
% voxels_baseline_mean = mean(voxels_data_4D_mask1(:,:,:,first_baseline_point_to_include:last_baseline_point_mask2),4);
% add "eps" to the denominator to avoid dividing 0:0 in the out-of-mask1
% voxels:
voxels_bl_min_ratio=abs(voxels_baseline_mean./(voxels_min_vals_in_time+eps));

th2=1.5; % th for relative value of (baseline average) / (min of bolus)
[row_data2,col_data2,slice_data2]=ind2sub(size(voxels_bl_min_ratio),find(voxels_bl_min_ratio>=th2));
[row_noise2,col_noise2,slice_noise2]=ind2sub(size(voxels_bl_min_ratio),find(voxels_bl_min_ratio<th2));
data2_indices=find(voxels_bl_min_ratio>=th2);
noise2_indices=find(voxels_bl_min_ratio<th2);
N_data2_voxels=length(data2_indices)

mask2=mask1;
mask2(noise2_indices)=0;

% applying Mask2 for the full 4D data:
voxels_data_4D_mask2=zeros(N_rows,N_cols,N_slices,N_time_points);
for t=1:N_time_points
   voxels_data_4D_mask2(:,:,:,t)=voxels_data_4D_mask1(:,:,:,t).*mask2;
end

% statistics for the data voxels after applying Mask2:
mean_DSC_curve=zeros(1,N_time_points);
std_DSC_curve=zeros(1,N_time_points);
for t=1:N_time_points
    voxels_data_3D_mask2_t=voxels_data_4D_mask2(:,:,:,t);
    voxels_data_3D_mask2_t_vec=voxels_data_3D_mask2_t(find(voxels_data_3D_mask2_t>0));
    mean_DSC_curve(t)=mean(voxels_data_3D_mask2_t_vec);
    std_DSC_curve(t)=std(voxels_data_3D_mask2_t_vec);
%     mean_DSC_curve(t)=sum(sum(sum(voxels_data_4D_mask1(:,:,:,t))))/N_data1_voxels;
%     mean_DSC_curve(t)=mean(voxels_data_4D(data1_indices));
end
%plot some random DSC curves and the mean of ALL DSC curves
N_DSC_curves_to_show=30;
check_vec=sort(random('unid',N_data2_voxels,1,N_DSC_curves_to_show));

if 1
    figure;
    for ii=1:length(check_vec)
        check_ind=check_vec(ii);
        plot(squeeze(voxels_data_4D_mask2(row_data2(check_ind),col_data2(check_ind),slice_data2(check_ind),:)));
        hold all;
    end
    std_part_to_show=1;
    title([num2str(N_DSC_curves_to_show),' random DSC curves, of voxels passed mask2. blue = mean of passed voxels + ',num2str(std_part_to_show),'std. green = mean of ALL voxel before Mask2']);
%     errorbar(mean_DSC_curve,std_part_to_show*std_DSC_curve,'--b','LineWidth',2);
    plot(mean_DSC_curve,'--b','LineWidth',3);
%     plot(mean_DSC_curve_total,'-.g','LineWidth',3);
    hold off;
end
%% statistics for the data after mask2 before mask3:

% time derivative / slope:
for t=1:N_time_points-1
   voxels_3D_diff_in_time(:,:,:,t)=voxels_data_4D_mask1(:,:,:,t+1)-voxels_data_4D_mask1(:,:,:,t);
   voxels_3D_time_diff_t_vec=reshape(voxels_3D_diff_in_time(:,:,:,t),1,N_rows*N_cols*N_slices);
   mean_DSC_time_diff_curve(t)=mean(voxels_3D_time_diff_t_vec);
   std_DSC_curve(t)=std(voxels_3D_time_diff_t_vec);
end


N_time_diff_curves_to_show=30;
check_vec=sort(random('unid',N_data1_voxels,1,N_time_diff_curves_to_show));

cc=hsv(N_time_diff_curves_to_show);
cc=cc(randperm(N_time_diff_curves_to_show),:);

% plot both time-diff curves and time curves:
if 1
   figure;
   std_part_to_show=1;
   for ii=1:length(check_vec)
        check_ind=check_vec(ii);
        plot(squeeze(voxels_data_4D_mask1(row_data1(check_ind),col_data1(check_ind),slice_data1(check_ind),:)),'color',cc(ii,:));
        hold on;
        plot(squeeze(voxels_3D_diff_in_time(row_data1(check_ind),col_data1(check_ind),slice_data1(check_ind),:)),'--','color',cc(ii,:));
        legendInfo{2*ii-1} = ['time curve #' num2str(ii)];
        legendInfo{2*ii} = ['time diff curve #' num2str(ii)];
        hold on;
   end
   plot(mean_DSC_curve,'--b','LineWidth',2);
   hold on;
   plot(mean_DSC_time_diff_curve,'--g','LineWidth',2);
   title([num2str(N_time_diff_curves_to_show),' random time-diff curves of DSC, of voxels passed mask1. blue = mean of passed voxels + ',num2str(std_part_to_show),'std. green = mean of time-difference ']);
   legendInfo{2*N_time_diff_curves_to_show+1}=['mean DSC time curve'];
   legendInfo{2*N_time_diff_curves_to_show+2}=['mean DSC time diff curve'];
   legend(legendInfo);
    
end


%% Attempts for threshold masking
im=voxels_data_4D(:,:,4,15);
% figure;
% for ii=0:11
%     th=ii*1000;
%     mask=zeros(size(im,1),size(im,2));
%     mask=double(im>=th);
% %     subplot(3,4,ii+1);
% %     imshow(mask);
% %     title(['ii=',num2str(ii)]);
%     im_th=im.*mask;
%     subplot(3,4,ii+1);
%     imshow(mat2gray(im_th));
%     title(['ii=',num2str(ii)]);
% end

