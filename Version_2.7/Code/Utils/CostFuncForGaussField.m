function Out=CostFuncForGaussField(x,pixels,SolVec)
mu=x(1:3);
L=zeros(3);
L([1 4 5 7 8 9])=x(4:9);
R= (L')*L;
invR=pinv(R);

Y1=pixels-ones(size(pixels,1),1)*mu;
Y2=-0.5*Y1*invR;
Y3 = dot(Y1,Y2,2);
Out=gCost(SolVec,Y3.^x(10),{'SumAbs' 'AddCompensate' 'MultiCompensate'});
% Out=gCost(SolVec,Y3,{'SumAbs' 'AddCompensate' 'MultiCompensate'});