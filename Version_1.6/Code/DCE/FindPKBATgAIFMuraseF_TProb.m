function [PKs Sims Cost] = FindPKBATgAIFMuraseF(DataToFit,SAIF,SampleTs,CumSumAIFdTSamp,WPerCTC)
% Find Kep using diffrence of C, as in Murase 2004
WarningS=warning('off','MATLAB:rankDeficientMatrix');
dT=SampleTs(2)-SampleTs(1);
% BAT Kep Vp Ktrans Ve RMS fVal f0Val
nToFit=size(DataToFit,1);
nTDif=size(SAIF,1);
nVols=size(DataToFit,2);
p0=0;
p1=2;
p2=4;
alpha=0.01;
FThresh = finv(1-alpha,p2-p1,nVols-p2);
F0Thresh = finv(1-alpha,p2-p0,nVols-p2);

PKs=zeros(nToFit,6);
Sims=zeros(size(DataToFit));

MinVp=0;
MaxVp=1;
MinVe=0.05;
MaxVe=10; % Account for AIF maxAmp too big -> ktrans small
MaxKtrans=2.2;
MaxKep=MaxKtrans/MinVe;
% MaxVes=min(1,2.2./Kep);
% KtransByT1 VeByT1 VpByT1 0 KepByT1

CForCluster=zeros(numel(SampleTs),nTDif);

VeByTdif=zeros(1,nTDif);
KtransByTdif=zeros(1,nTDif);
KepByTdif=zeros(1,nTDif);
VpByTdif=zeros(1,nTDif);
VpByTdifNoP=zeros(1,nTDif);

% mCumSumDatadT=cumtrapz(DataToFit,2).*dT;
mCumSumDatadT=-cumtrapz(SampleTs,DataToFit')';
% mCumSumDatadT=-[zeros(nToFit,1) cumsum(DataToFit(:,1:end-1),2)]; % *dT;
MeanAbsWDiff=NaN(1,nTDif);
MeanAbsWDiffNoP=NaN(1,nTDif);
CForCluster=NaN(nTDif,size(DataToFit,2));
% C for each T can be computed out of AIF Search in prepare

for c=1:nToFit
    for t=1:nTDif
        A=[CumSumAIFdTSamp(t,:)' mCumSumDatadT(c,:)' SAIF(t,:)'];
%         figure;plot(SampleTs,A(:,1),'k',SampleTs,A(:,2),'r')
%         figure;plot(SampleTs,A(:,1),'k',SampleTs,A(:,2),'b',SampleTs,DataToFit(c,:),'r',SampleTs,SAIF(t,:),'g')
        B=A\(DataToFit(c,:)'); % pinv(A)*CurC; % Find easy matrix inverse A\C
        if(B(2)<0)
            Tmp=A(:,[1 3])\(DataToFit(c,:)');
            B=[Tmp(1); 0; Tmp(2)];
        end
        if(B(1)-B(2)*B(3)<0)
            B=[0 0 0]';
        end
        % B is [Ktrans+Kep*Vp, Kep, Vp]
        KepByTdif(t)=min(MaxKep,max(eps,B(2)));
        VpByTdif(t)=min(MaxVp,max(MinVp,B(3)));
        KtransByTdif(t)=min(MaxKtrans,max(0,B(1)-B(2)*B(3)));
        VeByTdif(t)=max(MinVe,min(MaxVe,KtransByTdif(t)/KepByTdif(t)));
        KepByTdif(t)=KtransByTdif(t)/VeByTdif(t);
        
        CForCluster(t,:)=A*[B(1); KepByTdif(t); VpByTdif(t)];
        Diff=CForCluster(t,:)-DataToFit(c,:);
        MeanAbsWDiff(t)=sum((Diff.^2));
        
        % No permeability
        VpByTdifNoP(t)=min(MaxVp,max(MinVp,(SAIF(t,:)')\(DataToFit(c,:)')));
        CForClusterNoP(t,:)=(SAIF(t,:)')*VpByTdifNoP(t);
        Diff=CForClusterNoP(t,:)-DataToFit(c,:);
        MeanAbsWDiffNoP(t)=sum((Diff.^2));
        
%         figure;plot(SampleTs,DataToFit(c,:),'k',SampleTs,CForCluster(t,:),'b');
%         figure;plot(SampleTs,[DataToFit(c,:); Ax(:,1)'*B(1); Ax(:,2)'*B(2); Ax(:,3)'*B(3);CForCluster(t,:)])
    end
    % BAT Kep Vp Ktrans Ve RMS
    RSSNoModel=sum(DataToFit(c,:).^2);
    RSSNoP=min(MeanAbsWDiffNoP);
    RSSa=min(MeanAbsWDiff);
    fval=((RSSNoP-RSSa)/(p2-p1))/(RSSa/(nVols-p2));
    f0val=((RSSNoModel-RSSa)/(p2-p0))/(RSSa/(nVols-p2));
    PKs(c,7)=fval;
    PKs(c,8)=f0val;
    if(f0val<F0Thresh)
        PKs(c,6)=1;
        PKs(c,4)=0;
        PKs(c,5)=0;
        PKs(c,3)=0;
        PKs(c,2)=0;
        PKs(c,1)=1;
        Sims(c,:)=CForClusterNoP(1,:)*0;
        PKs(c,9)=2;
    else
        if(fval>FThresh)
            [PKs(c,6), MinIdxs]=min(MeanAbsWDiff);
            PKs(c,4)=KtransByTdif(MinIdxs);
            PKs(c,5)=VeByTdif(MinIdxs);
            PKs(c,3)=VpByTdif(MinIdxs);
            PKs(c,2)=KepByTdif(MinIdxs);
            PKs(c,1)=MinIdxs;
            Sims(c,:)=CForCluster(MinIdxs,:);
            
            Others=setdiff(1:nTDif,(MinIdxs-2:MinIdxs+2));
            if(~isempty(Others))
                PKs(c,9)=PKs(c,6)./min(MeanAbsWDiff(Others));
            end
        else
            [PKs(c,6), MinIdxs]=min(MeanAbsWDiffNoP);
            PKs(c,4)=0;
            PKs(c,5)=0;
            PKs(c,3)=VpByTdifNoP(MinIdxs);
            PKs(c,2)=0;
            PKs(c,1)=MinIdxs;
            Sims(c,:)=CForClusterNoP(MinIdxs,:);
            
            Others=setdiff(1:nTDif,(MinIdxs-2:MinIdxs+2));
            if(~isempty(Others))
                PKs(c,9)=PKs(c,6)./min(MeanAbsWDiffNoP(Others));
            end
        end
    end
end
warning(WarningS);
PKs(:,7)=1-fcdf(PKs(:,7),p2-p1,nVols-p2);
PKs(:,8)=1-fcdf(PKs(:,8),p2-p0,nVols-p2);

if(nargout>=3)
    Cost=max(sqrt(mean(((Sims-DataToFit).^2),2)).*WPerCTC);
end

if(false)
    figure(111);clf;
    nPlots=nToFit;
    for i=1:nPlots
        IStart=0;
        gsubplot(nPlots,i);
        plot(SampleTs,DataToFit(i+IStart,:),'k','LineWidth',2);hold on;plot(SampleTs, Sims(i+IStart,:),'b');
        title([num2str(PKs(i+IStart,1))]);
    end
end