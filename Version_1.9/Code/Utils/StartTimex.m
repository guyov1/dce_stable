img1=loadniidata([WorkingP 'DCE4D.nii']);
img2=img1(:,:,:,:);
tmp=mean(img2(:,:,:,1:4),4);
s=size(img2);
normMat=zeros(s(1),s(2),s(3),s(4));

for i=1:s(4)
normMat(:,:,:,i)=tmp;
end;
normMat=double(img2)*100./double(normMat)-100;
[A B C IMap]=ndgrid(1:s(1),1:s(2),1:s(3),1:s(4));

%%
% threshMat=normMat>10;
% BnormMat=BnormMat;
tmp=(normMat>5).*IMap;
B=tmp*0;
B(tmp==0)=1;
B(:,:,:,1:4)=1;
tmp1=tmp;
tmp(B>0)=s(4);
startTime=min(tmp,[],4);
startTime(max(tmp1,[],4)==0)=0;

% startTime=zeros(s(1),s(2),s(3));
% for i=1:s(1)
%     for j=1:s(2)
%         for k=1:s(3)
%             vecTemp=squeeze(normMat(i,j,k,:));
%             ind=find(vecTemp>10,1,'first');
%             if length(ind)>0
%                 startTime(i,j,k)=ind;
%             end;
%                 
%         end;
%     end;
% end;
mricront(startTime)
% Raw2Nii(startTime,B1FN,'float32', MeanFN);
%%
% CTCFN=[WorkingP 'AfterCTC' USStr '.mat'];
% load(CTCFN,'CTC2DBigGood','MskCTCGood');
ACTC2D=CTC2DBigGood(:,1:20);
ACTC2D(:,1:4)=0;
NCTC2D=ACTC2D./repmat(max(ACTC2D,[],2),[1 size(ACTC2D,2)]);
[A, IMapx]=ndgrid(1:size(ACTC2D,1), 1:size(ACTC2D,2));
tmp=100-99*(NCTC2D>0.1);
startTime2D=min(tmp.*IMapx,[],2);
startTime2D(startTime2D==100)=0;

startTime=MskCTCGood3D*0;
startTime(MskCTCGood3D)=startTime2D;
mricront(startTime)

