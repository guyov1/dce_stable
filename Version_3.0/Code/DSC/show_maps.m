% read and show maps

cd maps;
% Flow:
load('Flow_sSVD.mat');
F_sSVD=F;
N_slices=size(F_sSVD,3);

figure;
n_sp_rows=floor(sqrt(N_slices));
n_sp_cols=ceil(N_slices/n_sp_rows);
for slice=1:N_slices
    subplot(n_sp_rows,n_sp_cols,slice);
    imshow(mat2gray(F_sSVD(:,:,slice)));
end