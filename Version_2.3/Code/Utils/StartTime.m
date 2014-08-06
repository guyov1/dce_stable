data=load_Untouch_nii('dsc4drm.nii');
img1=data.img;
img2=img1(:,:,:,:);
tmp=mean(img2(:,:,:,1:6),4);
s=size(img2);
normMat=zeros(s(1),s(2),s(3),s(4));

for i=1:s(4)
normMat(:,:,:,i)=tmp;
end;
normMat=double(img2)*100./double(normMat)-100;


threshMat=normMat<-5;
startTime=zeros(s(1),s(2),s(3));
for i=1:s(1)
    for j=1:s(2)
        for k=1:s(3)
            vecTemp=squeeze(normMat(i,j,k,:));
            ind=find(vecTemp<-10,1,'first');
            if length(ind)>0
                startTime(i,j,k)=ind;
            end;
                
        end;
    end;
end;
data2=load_Untouch_nii('vol1.nii');
startTimeDSC=data2;
startTimeDSC.img=startTime;
save_untouch_nii(startTimeDSC,'startTimeDSC.nii');
