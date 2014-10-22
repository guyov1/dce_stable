function Prepared=findPK_1LevelPrepareKnownT1(ParamMat,TimePoints,RCTCE,TimeBetweenDCEVols,TWM,CWs,SampleFactor)

% ParamMat - [ min max power nPoints_in_1st_Level nPoints_in_other_Levels ]

ParamsV=cell(1,size(ParamMat,1));
for i=1:size(ParamMat,1)
    ParamsV{i}=gpowspace(ParamMat(i,1),ParamMat(i,2),ParamMat(i,4),ParamMat(i,3));
end

[ParamsM{1:3}]=ndgrid(ParamsV{:});
Params=[ParamsM{1}(:) ParamsM{2}(:) ParamsM{3}(:)];
Prepared.ParamsV=ParamsV; % for debug only
Prepared.Params=Params;
Prepared.ParamMat=ParamMat;
Prepared.SampleFactor=SampleFactor;

% NRCDCE=repMulti(RCDCE,1./max(RCDCE,[],2));
Prepared.BizzareCs=[]; % FindBizzare(mean(abs(diff(NRCDCE(:,(BolusStart+10):end),1,2)),2),1);
% disp(['Ignoring clusters ' num2str(Prepared.BizzareCs')]);
Prepared.NotBizzareCs=setdiff(1:size(RCTCE,1),Prepared.BizzareCs);
RCTCE=RCTCE(Prepared.NotBizzareCs,:);

nTimePointsForPKsearchExh=numel(TimePoints);
TimePointsForPKsearchExhMin=(TimePoints-1)*TimeBetweenDCEVols/60;
[X Y]=meshgrid(1:nTimePointsForPKsearchExh);
% ConvIdxM=(repmat(TimePoints,nTimePoints,1).*(X<=Y));
ConvIdxM=((Y-X)+1).*(X<=Y);

Prepared.nCClusters=size(RCTCE,1);
SamplePs=1:SampleFactor:nTimePointsForPKsearchExh;
OrigTimePointsMin=TimePointsForPKsearchExhMin(SamplePs);

nTimePointsOrig=numel(SamplePs);
dTOrig=OrigTimePointsMin(2)-OrigTimePointsMin(1);
[X Y]=meshgrid(1:nTimePointsOrig);
% ConvIdxM=(repmat(TimePoints,nTimePoints,1).*(X<=Y));
ConvIdxMOrig=((Y-X)+1).*(X<=Y);
Prepared.TriBOrig=ConvIdxMOrig>0;
Prepared.AllConvIdxMOrig=ConvIdxMOrig(Prepared.TriBOrig); % AllConvIdxMOrig;

% TW=OrigTimePointsMin*0;
% TW(1:BolusStart-1)=1/(BolusStart-1);
% TW((BolusStart-1):(BolusStart+2))=3/3;
% TW(BolusStart)=3;
% TW((BolusStart+2):end)=1/(nTimePointsForPKsearchExh-(BolusStart+1));
TW=OrigTimePointsMin*0+1;
% if(isempty(TWM))
%     Prepared.TWM=repmat(TW,[size(Params,1),1,Prepared.nCClusters]);
% else
%     Tmp=repMulti(TWM,TW);
%     Prepared.TWM=repmat(permute(Tmp,[3 2 1]),[size(Params,1) 1 1]);
% end

% Kep=Params(:,1)./Params(:,2);
% [U, ~, IB]=unique(Kep);
% nUKep=numel(U);
% Prepared.EKepM=cell(1,nUKep);
% TriB=ConvIdxM>0;
% for i=1:nUKep
%     EKepV=exp(-U(i)*TimePointsForPKsearchExhMin);
%     Prepared.EKepM{i}=zeros(size(ConvIdxM));
%     Prepared.EKepM{i}(TriB)=EKepV(ConvIdxM(TriB));
% end

% Prepared.IdxMat=repmat(IB-1,1,nTimePointsForPKsearchExh)*nTimePointsForPKsearchExh+repmat(1:nTimePointsForPKsearchExh,size(Params,1),1);
% E1firstLevel=exp(-TR./Params(:,4));
% E1MfirstLevel=repmat(E1firstLevel,1,nTimePointsForPKsearchExh);
% E1MfirstLevelCosFA=E1MfirstLevel*Prepared.CosFA;
% Prepared.B1M=(1-E1MfirstLevelCosFA)./(1-E1MfirstLevel);
% Prepared.R10M=repmat(1./Params(:,4),1,nTimePointsForPKsearchExh);
% Prepared.KtransM=repmat(Params(:,1),1,nTimePointsForPKsearchExh);
% Prepared.RCDCEC=cell(1,Prepared.nCClusters);
Prepared.RCTCE=RCTCE;
% for Pyr_i=1:Prepared.nCClusters
%     Prepared.RCDCEC{Pyr_i}=repmat(-RCDCE(Pyr_i,OrigTimePoints),size(Params,1),1);
% end
Prepared.AllTimePointsMinute=TimePointsForPKsearchExhMin; % AllTimePointsMinute;
Prepared.TriB=ConvIdxM>0;
Prepared.AllConvIdxM=ConvIdxM(Prepared.TriB); % AllConvIdxM;
Prepared.TW=TW;
if(numel(CWs)<=1)
    CWs=ones(1,Prepared.nCClusters);
end
Prepared.CWs=CWs./sum(CWs);
% Keps=unique(gpowspace(0,2.2,11,5)'*(1./gpowspace(0.15,1,11,4)));
% if(numel(Keps)>nOptsKepT1(1))
%     A=GaussianMixture(Keps,nOptsKepT1(1),nOptsKepT1(1),false);
%     KepsC=[A.cluster.mu];
% else
%     KepsC=Keps;
% end
% Prepared.KepsC=KepsC';

% %%
% T1s=200:10:5000;
% C=0:0.0001:0.04;
% RM=zeros(numel(C),700);
% for i=1:numel(C)
%     for j=1:numel(T1s)
%         CurR=RfuncC(1./T1s(j),C(i),TR,Prepared.CosFA);
%         CurRIdx=floor(CurR*10);
%         RM(i,CurRIdx)=T1s(j);
%     end
% end
% % Rs=(1:700)/10;
% % figure(75858);surf(Rs(1:10:end),C(1:10:end),RM(1:10:end,1:10:end))
% Prepared.T1ByCRmat=RM;

Prepared.OrigTimePointsMin=OrigTimePointsMin;
Prepared.SamplePs=SamplePs;
Prepared.nTimePointsOrig=nTimePointsOrig;
Prepared.dTOrig=dTOrig;

% for given T1
%Cs

% R10=1./ParamsV{4};
% nT1s=numel(R10);
% E1=exp(-TR*R10);
% Bt=repmat(RCDCE,[1 1 nT1s]).*repmat(permute((1-E1)./(1-E1*Prepared.CosFA),[1 3 2]),[size(RCDCE) 1]);
% 
% Et=(1-Bt)./(1-Bt*Prepared.CosFA);
% C=-log(Et)/TR-repmat(permute(R10,[1 3 2]),[size(RCDCE) 1]);
% Prepared.C=permute(C,[2 1 3]);
% 
% Prepared.R102D=repmat(R10,[nTimePointsOrig,1]);
% OneMinusE1CosFAdOneMinusE1=(1-E1*Prepared.CosFA)./(1-E1);
% Prepared.OneMinusE1CosFAdOneMinusE12D=repmat(OneMinusE1CosFAdOneMinusE1,[nTimePointsOrig,1]);
% Prepared.TWRep=repmat(TW',[1 nT1s]);
% Prepared.nT1s=nT1s;
% Prepared.R10=R10;
