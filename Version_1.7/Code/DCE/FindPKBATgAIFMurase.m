function [PKs Sims Cost] = FindPKBATgAIFMurase(WorkingP,DataToFit,SAIF,SampleTs,CumSumAIFdTSamp,WPerCTC)

% Find Kep using diffrence of C, as in Murase 2004
WarningS=warning('off','MATLAB:rankDeficientMatrix');

% Get the time derivative 
dT=SampleTs(2)-SampleTs(1);

%% BAT Kep Vp Ktrans Ve RMS calculation

% Number of points to fit
nToFit=size(DataToFit,1);
% Number of shifted AIFs
nTDif=size(SAIF,1);

% Initiate paremters
PKs=zeros(nToFit,6);
Sims=zeros(size(DataToFit));

% MinVp=0;
% MaxVp=1;
% MinVe=0.05;
% MaxVe=10;
% MaxKtrans=2.2;
% MaxKep=MaxKtrans/MinVe;

% Some parameters usefull for better fitting
MinVp=0;
MaxVp=1*10; % Account for AIF maxAmp too small
MinVe=0.05;
MaxVe=10; % Account for AIF maxAmp too big -> ktrans small
MaxKtrans=2.2*10; % Account for AIF maxAmp too small
MaxKep=MaxKtrans/MinVe;
% MaxVes=min(1,2.2./Kep);
% KtransByT1 VeByT1 VpByT1 0 KepByT1

CForCluster=zeros(numel(SampleTs),nTDif);

VeByTdif=zeros(1,nTDif);
KtransByTdif=zeros(1,nTDif);
KepByTdif=zeros(1,nTDif);
VpByTdif=zeros(1,nTDif);

% mCumSumDatadT=cumtrapz(DataToFit,2).*dT;

% Use the time vector for cumtrapz
CTC_FN=[WorkingP 'AfterCTC.mat'];

Num_Of_T1_Maps=0;
% load(CTC_FN,'Num_Of_T1_Maps'); % Guy: change this. Takes time, it's
% called a lot, it's inside fminsearch

% Make sure we actually have more than 1 T1 base map
if (Num_Of_T1_Maps>1)
    error('Change above!');
    load(CTC_FN,'Complete_Time_Vector');
    Complete_Time_Vector_Minutes = Complete_Time_Vector/60;

    mCumSumDatadT=-cumtrapz(Complete_Time_Vector_Minutes,DataToFit')';

else  % Regular case, 1 T1 map, no time vector

    mCumSumDatadT=-cumtrapz(SampleTs,DataToFit')';
end

% mCumSumDatadT=-[zeros(nToFit,1) cumsum(DataToFit(:,1:end-1),2)];
MeanAbsWDiff=NaN(1,nTDif);

CForCluster=NaN(nTDif,size(DataToFit,2));
% C for each T can be computed out of AIF Search in prepare

% Go over all needed voxels
for c=1:nToFit
     % For each possible shifted AIF
    for t=1:nTDif
        A=[CumSumAIFdTSamp(t,:)' mCumSumDatadT(c,:)' SAIF(t,:)'];
        B=A\(DataToFit(c,:)'); % pinv(A)*CurC; % Find easy matrix inverse A\C

        % If kep < 0, calculate again ignoring it
        if(B(2)<0)
            Tmp=A(:,[1 3])\(DataToFit(c,:)');
            B=[Tmp(1); 0; Tmp(2)];
        end
        
        % B is [Ktrans+Kep*Vp, Kep, Vp]
        KepByTdif(t)=min(MaxKep,max(0,B(2)));
        VpByTdif(t)=min(MaxVp,max(MinVp,B(3)));
        KtransByTdif(t)=min(MaxKtrans,max(0,B(1)-B(2)*B(3)));
        VeByTdif(t)=max(MinVe,min(MaxVe,KtransByTdif(t)/KepByTdif(t)));
        KepByTdif(t)=KtransByTdif(t)/VeByTdif(t);
        
        if(B(1)-B(2)*B(3)<0)
            B=(SAIF(t,:)')\(DataToFit(c,:)');
            KtransByTdif(t)=0;
            KepByTdif(t)=0;
            VeByTdif(t)=0;
            VpByTdif(t)=min(MaxVp,max(MinVp,B(1)));
            B(1)=0;
        end
        
%         Ax=[CumSumAIFdTSamp(t,:)' mCumSumDatadT(c,:)' SAIF(t,:)'];
        CForCluster(t,:)=A*[B(1); KepByTdif(t); VpByTdif(t)];
        Diff=CForCluster(t,:)-DataToFit(c,:);
        MeanAbsWDiff(t)=sqrt(sum((Diff.^2)'));
        
%         figure;plot(SampleTs,DataToFit(c,:),'k',SampleTs,CForCluster(t,:),'b');
%         figure;plot(SampleTs,[DataToFit(c,:); Ax(:,1)'*B(1); Ax(:,2)'*B(2); Ax(:,3)'*B(3);CForCluster(t,:)])
    end
    % BAT Kep Vp Ktrans Ve RMS
    [PKs(c,6), MinIdxs]=min(MeanAbsWDiff);
    PKs(c,4)=KtransByTdif(MinIdxs);
    PKs(c,5)=VeByTdif(MinIdxs);
    PKs(c,3)=VpByTdif(MinIdxs);
    PKs(c,2)=KepByTdif(MinIdxs);
    PKs(c,1)=MinIdxs;
    
    Sims(c,:)=CForCluster(MinIdxs,:);
end
warning(WarningS);

if(nargout>=3)
    Diff=Sims-DataToFit;
    SqrDiff=(Diff).^2;
%     AllCosts=sqrt(mean(SqrDiff,2));
    AllCosts=sqrt(max(SqrDiff,[],2));
    Cost=mean(AllCosts.*WPerCTC);
%     Cost=max(AllCosts.*WPerCTC);
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
a=5;