clear all;
%% Dicom source path:
%SRC_DIR='D:\users\chenh\DSC_project\files_from_guy\DCE healthy data\sub10APOL_SHIRA\DCE\St09_Se16_DCE with FA 25';
%TRG_FN='D:\users\chenh\DSC_project\temp\output5';
%% Dicom to NiFTy:
%gDicom2Nifti(SRC_DIR,TRG_FN);
%%
nii_path='D:\users\chenh\DSC_project\temp\output5\ArSi_20081130\DCEMainCoreged\Coreged_vol_0020.nii';
[Out Header A]=loadniidata(nii_path);

figure;
im=diff;
num_images=size(im,3);
num_images_per_row=5;
for ii=3:3
    image_ii=mat2gray(im(:,:,ii));
%     subplot(ceil(num_images/num_images_per_row),num_images_per_row,ii);
    imshow(image_ii);
end

figure;
im=A.img;
num_images=size(im,3);
num_images_per_row=5;
for ii=3:3
    image_ii=mat2gray(im(:,:,ii));
%     subplot(ceil(num_images/num_images_per_row),num_images_per_row,ii);
    imshow(image_ii);
end