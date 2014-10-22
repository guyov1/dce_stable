%%
a=GaussianMixture([I{1} I{2} I{3}],1,1,false);
L=chol(a.cluster.R);
x0=[a.cluster.mu' L([1 4 5 7 8 9]) 1];
pixels=[I{1} I{2} I{3}];
CostFunc=@(x) CostFuncForGaussField(x,pixels,SolVec);
%     BestX=fminsearch(CostFunc,x0(1:9));
BestX=fminsearch(CostFunc,x0);
x=BestX;
%%
mu=x(1:3);
L=zeros(3);
L([1 4 5 7 8 9])=x(4:9);
R= (L')*L;
invR=pinv(R);

Y1=pixels-ones(size(pixels,1),1)*mu;
Y2=-0.5*Y1*invR;
Y3 = dot(Y1,Y2,2);
% Out=gCost(SolVec,Y3.^x(10),'SumAbs',{'AddCompensate' 'MultiCompensate'});
%     [Out More k]=gCost(SolVec,Y3,{'SumAbs' 'AddCompensate' 'MultiCompensate'});
[Out More k]=gCost(SolVec,Y3.^x(10),{'SumAbs' 'AddCompensate' 'MultiCompensate'});

Y1F=[IF{1} IF{2} IF{3}]-ones(size([IF{1} IF{2} IF{3}],1),1)*mu;
Y2F=-0.5*Y1F*invR;
Y3F = dot(Y1F,Y2F,2).^x(10);
NewField(:)=abs(Y3F*k(1)+k(2));
CT1=(UncleanedT1./NewField)*RefTrgVal;
sqrt(mean(((CT1(RefVolX)-RefTrgVal)).^2))
sqrt(mean(W.*((CT1(RefVolX)-RefTrgVal)).^2))