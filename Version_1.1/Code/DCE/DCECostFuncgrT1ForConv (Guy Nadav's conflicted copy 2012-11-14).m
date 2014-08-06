% function Out=DCECostFunc(AIF,Ktrans,Ve,Vp,TimePoints,ConvIdxMTriB,TR,CosFA,R10,r1,Ws,RData)
function Out=DCECostFuncgrT1ForConv(AIF,Kep,TimePoints,ConvIdxMTriB,TriB)
n=numel(Kep);
nTimePoints=numel(TimePoints);
% CosFA=cosd(FA);

dT=TimePoints(2)-TimePoints(1);
[U, QQQ, IB]=unique(Kep);
nUKep=numel(U);
EKepM=zeros(nTimePoints);
% TriB=ConvIdxM>0;
CAIF=zeros(nTimePoints,nUKep);
for i=1:nUKep
    EKepV=exp(-U(i)*TimePoints);
%     EKepM(TriB)=EKepV(ConvIdxM(TriB));
    EKepM(TriB)=EKepV(ConvIdxMTriB);
    CAIF(:,i)=EKepM*AIF;
end
IdxMat=repmat(IB-1,1,nTimePoints)*nTimePoints+repmat(1:nTimePoints,n,1);
% get C
if(numel(IB)==1)
    Out=CAIF(IdxMat)'*dT;
else
    Out=CAIF(IdxMat)*dT;
end