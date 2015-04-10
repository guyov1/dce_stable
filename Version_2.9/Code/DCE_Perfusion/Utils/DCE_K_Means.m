%clear all
File_4D_Path    = 'D:\users\guyn\DCE_OUT\Guy_Test_1\SmVl_20120930\DCE4D.nii';
Brain_Mask_Path = 'D:\users\guyn\DCE_OUT\Guy_Test_1\SmVl_20120930\BrainMask.nii';

DCE_4D      = loadniidata(File_4D_Path); 
Brain_Mask  = loadniidata(Brain_Mask_Path);



blabla           = squeeze(DCE_4D(:,:, :,:));

[m n d s ]= size(blabla);

N = m*n*d ;

rsnormMat               = reshape(blabla, N, s);
rsnormMat(rsnormMat==0) = NaN;

[IDX, C, sumd, D]       = kmeans(rsnormMat,5, 'Distance', 'cosine','EmptyAction','error','Start','cluster','Replicates',5);

m = 256;
n = 256;

[v1 s1]=size(rsnormMat);

d1     = v1 / (m*n);
IDXmat = reshape(IDX(1:v1), m, n, d1);


KmeansMap     = datax;
%KmeansMap.img = IDXmat;
    
save_untouch_nii(KmeansMap,'DCEKMout.nii');

