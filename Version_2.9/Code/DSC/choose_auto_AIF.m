function AIF=choose_auto_AIF(concentration_4D,Mask)

[row,col,slice]=ind2sub(size(concentration_4D),find(Mask>0));
N_data_voxels=length(row);
%% choose voxels randomally: (in the future it should be changed)
N_voxels_AIF=10;
% weight vector for the N_voxels_AIF concentration curves
W = (1/N_voxels_AIF) * ones(1,N_voxels_AIF);
voxels_AIF_vec=sort(random('unid',N_data_voxels,1,N_voxels_AIF));
AIF_voxels_curves=zeros(N_voxels_AIF,size(concentration_4D,4));
figure;hold all;
for ii=1:N_voxels_AIF
    voxel_AIF_ind=voxels_AIF_vec(ii);
    AIF_voxels_curves(ii,:)=concentration_4D(row(voxel_AIF_ind),col(voxel_AIF_ind),slice(voxel_AIF_ind),:);
    plot(AIF_voxels_curves(ii,:));
end
% the AIF is a weighted sum of some voxels concentration curves:
AIF=W*AIF_voxels_curves;
figure;plot(AIF);title(['auto AIF chosen, from ',num2str(N_voxels_AIF),' random voxels that passed pre-processing filters'],'FontSize',13); 
%%

end