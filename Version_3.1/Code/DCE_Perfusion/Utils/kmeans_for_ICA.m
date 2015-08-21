function [] = kmeans_for_ICA(NiiPath, NiiName)


%NiiPath  = 'DCEmainBrain.nii';
NiiImg   = [NiiPath filesep NiiName];
data     = load_Untouch_nii(NiiImg);
dataimg  = data.img;
InVolume = squeeze(dataimg(:,:, :,:));

%size(x)
%imagesc(x(:,:,8))
%plot(squeeze(data(130,90,4,:)));%x=90;y=130

tmp      = mean(InVolume(:,:,:,1:3),4);
s        = size(InVolume);
normMat  = zeros(s(1),s(2),s(3),s(4));

for i=1 : s(4)
    normMat(:,:,:,i)=tmp;
end;

xnormMat = InVolume*100./normMat-100;
imagesc(xnormMat(:,:,4,11))
%plot(squeeze(normMat(130,90,11,:)));

[m, n, d, s]      = size(xnormMat);
N                 = m*n*d;
rsnormMat         = reshape(xnormMat, N, s);
[IDX, C, sumd, D] = kmeans(rsnormMat,5);

IDXmat            = reshape(IDX, m, n, d);
imagesc(IDXmat(:,:,5))

figure
plot(C');
legend('c1','c2','c3','c4','c5','c6','c7')

datax         = load_Untouch_nii('ref.nii');
KmeansMap     = datax;
KmeansMap.img = IDXmat;
Output_Nii    = [NiiPath filesep 'KmeansMap.nii'];
save_untouch_nii(KmeansMap,Output_Nii);

end



