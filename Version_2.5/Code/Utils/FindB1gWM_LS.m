% PolyDeg=2;
% Tmp=dec2base((0:(PolyDeg+1)^3-1),PolyDeg+1)-'0';
% Powers=[0 1 2 3 4];
% GoodOnes=Tmp((sum(Tmp,2)<=PolyDeg),:);
% GoodOnes=Powers(GoodOnes+1);
% Tmp=double(RefVolX);
% for i=1:size(RefVolX,3)
%     H = fspecial('gaussian',10,2);
%     RefVolXs(:,:,i)=imfilter(squeeze(Tmp(:,:,i)),H,'replicate');
% end
% clear I IF ProbMat FieldMat
% [I{1:3}]=ind2sub(size(RefVolX),find(RefVolX));
% [IF{1:3}]=ind2sub(size(RefVolX),find(RefVolX+1));
% % NonLinTransform=@(x) x.^a;
% % rNonLinTransform=@(x) x.^(1/a);
% %     for i=1:3
% %         I{i}=I{i}-size(RefVolX,i)/2;
% %         IF{i}=IF{i}-size(RefVolX,i)/2;
% %     end
% for i=1:size(GoodOnes,1)
%     ProbMat(:,i)=(I{1}.^GoodOnes(i,1)).*(I{2}.^GoodOnes(i,2)).*(I{3}.^GoodOnes(i,3));
%     FieldMat(:,i)=(IF{1}.^GoodOnes(i,1)).*(IF{2}.^GoodOnes(i,2)).*(IF{3}.^GoodOnes(i,3));
% end
% %     ProbMat=[I.^2 J.^2 K.^2 I.*J I.*K J.*K I J K I*0+1];
% % Func=@(x) [x(:,1).^0 x(:,1) x(:,2) x(:,3) x(:,1).*x(:,2) x(:,1).*x(:,3) x(:,2).*x(:,3)];
% % ProbMat=Func([I{1} I{2} I{3}]);
% % FieldMat=Func([IF{1} IF{2} IF{3}]);
% NonLinTransformPower=1;
% NonLinTransform=@(x) log(x.^(1/NonLinTransformPower));
% rNonLinTransform=@(x) exp(x).^NonLinTransformPower;
% UncleanedData=UncleanedT1(RefVolX);
% SolVec=NonLinTransform(UncleanedData);
% 
% W=(1./RefVolXs(RefVolX)).^1;
% 
% CostFunc=@(x) FindB1gWMCoeffsWithNonLinTransform(abs(x),UncleanedData,ProbMat,W);
% 
% [Cost ParamVec]=CostFunc(1);
% 
% 
% NewFieldVec=rNonLinTransform(FieldMat*ParamVec);
% NewField=RefVolX*0;
% NewField(:)=NewFieldVec/RefTrgVal;

% NonLinTransformPower=1;
% NonLinTransform=@(x) log(x.^(1/NonLinTransformPower));
% rNonLinTransform=@(x) exp(x).^NonLinTransformPower;
[NewField W I IF SolVec]=FitNonLinPolyField(2,NonLinTransform,rNonLinTransform,UncleanedT1,RefVolX);
NewField=NewField/RefTrgVal;

PseudoB1=sqrt(NewField);
CT1=(UncleanedT1./NewField);
B1gWMrRMS=sqrt(mean(W.*((CT1(RefVolX)./RefTrgVal-1)).^2));

AddToLog(WorkingP,'c_30B1gWMrRMS',['$B_1$ given WM rRMS: ' num2str(B1gWMrRMS)]);