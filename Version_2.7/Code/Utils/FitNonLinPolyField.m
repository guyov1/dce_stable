function [Field W I IF SolVec]=FitNonLinPolyField(PolyDeg,NonLinTransform,rNonLinTransform,DataToFit,Mask)
% PolyDeg=2;
Tmp=dec2base((0:(PolyDeg+1)^3-1),PolyDeg+1)-'0';
Powers=[0 1 2 3 4];
GoodOnes=Tmp((sum(Tmp,2)<=PolyDeg),:);
GoodOnes=Powers(GoodOnes+1);
Tmp=double(Mask);
for i=1:size(Mask,3)
    H = fspecial('gaussian',10,2);
    Masks(:,:,i)=imfilter(squeeze(Tmp(:,:,i)),H,'replicate');
end
clear I IF ProbMat FieldMat
[I{1:3}]=ind2sub(size(Mask),find(Mask));
[IF{1:3}]=ind2sub(size(Mask),find(Mask+1));
%     for i=1:3
%         I{i}=I{i}-size(Mask,i)/2;
%         IF{i}=IF{i}-size(Mask,i)/2;
%     end
for i=1:size(GoodOnes,1)
    ProbMat(:,i)=(I{1}.^GoodOnes(i,1)).*(I{2}.^GoodOnes(i,2)).*(I{3}.^GoodOnes(i,3));
    FieldMat(:,i)=(IF{1}.^GoodOnes(i,1)).*(IF{2}.^GoodOnes(i,2)).*(IF{3}.^GoodOnes(i,3));
end
%     ProbMat=[I.^2 J.^2 K.^2 I.*J I.*K J.*K I J K I*0+1];
% Func=@(x) [x(:,1).^0 x(:,1) x(:,2) x(:,3) x(:,1).*x(:,2) x(:,1).*x(:,3) x(:,2).*x(:,3)];
% ProbMat=Func([I{1} I{2} I{3}]);
% FieldMat=Func([IF{1} IF{2} IF{3}]);
UncleanedData=DataToFit(Mask);
SolVec=NonLinTransform(UncleanedData);

W=(1./Masks(Mask)).^1;

% CostFunc=@(x) FindB1gWMCoeffsWithNonLinTransform(UncleanedData,ProbMat,W,NonLinTransform,rNonLinTransform);

[Cost ParamVec]=FindB1gWMCoeffsWithNonLinTransform(UncleanedData,ProbMat,W,NonLinTransform,rNonLinTransform);

FieldVec=rNonLinTransform(FieldMat*ParamVec);
Field=Mask*0;
Field(:)=FieldVec;