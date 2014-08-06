function Msk3D=ShowImageWithPoints(Handles,Image,Msk,Idx,GoodRows,GoodCols,GoodSlices)

Msk3D=Msk*0>1;
Tmp2=find(Msk)*0>1;
Tmp2(Idx)=true;
Msk3D(Msk)=Tmp2;

TmpD=imdilate(Msk3D,strel('disk',4));
for i=1:numel(GoodSlices)
    I=squeeze(Image(:,:,GoodSlices(i)));
    I=min(4000,I)/4000;
    Tmp2=squeeze(TmpD(:,:,GoodSlices(i)));
    I(Tmp2)=0;
    Tmp2=squeeze(Msk3D(:,:,GoodSlices(i)));
    I(Tmp2)=1;
    IRGB(:,:,1)=I;
    I(Tmp2)=0;
    IRGB(:,:,2)=I;
    Tmp2=squeeze(TmpD(:,:,GoodSlices(i)));
    I(Tmp2)=1;
    Tmp2=squeeze(Msk3D(:,:,GoodSlices(i)));
    I(Tmp2)=0;
    IRGB(:,:,3)=I;
    IRGB3(:,:,:,i)=IRGB(GoodRows,GoodCols,:);
end
gfig(Handles);clf;
montage(mritransform(IRGB3))