function IRGB3=ShowImageWithPointsClr(Image,Msk,Idx,GoodRows,GoodCols,GoodSlices,ClrVec)

ClrM=[0 0 0; 1 0 0; 0 0 1; 0 1 0; 1 1 0; 1 0 1; 0 1 1];
DAll=double(Msk*0);
clear Msk3D
for i=1:max(ClrVec)
    Tmp2=find(Msk)*0>1;
    Tmp2(Idx(ClrVec==i))=true;
    Msk3D{i}=Msk*0>1;
    Msk3D{i}(Msk)=Tmp2;
    TmpD{i}=imdilate(Msk3D{i},strel('disk',4))*i;
    DAll(TmpD{i}~=0)=TmpD{i}(TmpD{i}~=0);
end
Tmp2=find(Msk)*0>1;
Tmp2(Idx)=true;
Msk3D=Msk*0>1;
Msk3D(Msk)=Tmp2;

for i=1:numel(GoodSlices)
    I=squeeze(Image(:,:,GoodSlices(i)));
    I=min(4000,I)/4000;
    IRGB=repmat(I,[1 1 3]);
    DRGB=ind2rgb(squeeze(DAll(:,:,GoodSlices(i)))+1,ClrM);
    BArea=repmat(any(DRGB~=0,3),[1 1 3]);
    IRGB(BArea)=DRGB(BArea);
    Tmp3=repmat(squeeze(Msk3D(:,:,GoodSlices(i))),[1 1 3]);
    BArea=repmat(any(Tmp3~=0,3),[1 1 3]);
    IRGB(BArea)=Tmp3(BArea);

    IRGB3(:,:,:,i)=IRGB(GoodRows,GoodCols,:);
end
IRGB3=mritransform(IRGB3);