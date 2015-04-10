P1='D:\users\Moran\OptDCE\BrainT\04_GOLDSHTEIN_CHEN_BELA\Study20140723_133944_T2\DCE\long\GoBe_20140723\';
P2='D:\users\Moran\OptDCE\BrainT\04_GOLDSHTEIN_CHEN_BELA\Study20140723_133944_T2\DCE\long - Copy\GoBe_20140723\';
art3d=loadniidata([P1 'InspectedRepVox.nii']);
CTC4D1=loadniidata([P1 'CTC4D.nii']);
CTC4D2=loadniidata([P2 'CTC4D.nii']);
art2D1=Reshape4d22d(CTC4D1,art3d);
art2D2=Reshape4d22d(CTC4D2,art3d);
PKMFN1=[P1 'PKM' '.mat'];
aa=load(PKMFN1);
load([P1 
figure;subplot(1,2,1);plot(aa.GoodTs,art2D1');
subplot(1,2,2);plot(aa.GoodTs,art2D2');