function [MIdxs CXs Sims] = FindKepBATgAIF(DataToFit,SAIF,SHConvd,Keps1I)
% BAT Kep Vp Ktrans
nToFit=size(DataToFit,1);
nKeps1=size(SHConvd,2);
nTDif=size(SAIF,1);
RMSs=zeros([nTDif nKeps1 nToFit]);
Xs=zeros([nTDif nKeps1 2 nToFit]);
WarningS=warning('off','MATLAB:rankDeficientMatrix');
for i=1:nTDif
    %     disp([i nTDif]);
    Regressors=[SAIF(i,:); SAIF(i,:)];
    for j=1:nKeps1
        Regressors(2,:)=SHConvd(i,Keps1I(j),:);
        % DataToFit=X*Regressors
        X=Regressors'\DataToFit';

        Xs(i,j,:,:)=X;
        Difs=((Regressors')*X)'-DataToFit;
        RMSs(i,j,:)=sqrt(mean(Difs.^2,2));
    end
end
CXs=zeros([2 nToFit ]);
MIdxs=zeros(nToFit,2);
for c=1:nToFit
    [Tmp MIdxs(c,:)]=gmin(squeeze(RMSs(:,:,c)));
    CXs(:,c)=squeeze(Xs(MIdxs(c,1),MIdxs(c,2),:,c));
end
FNeg=find(any(CXs<0,1));
v=ver('MATLAB');
if(str2num(v.Release(5:6))>10)
    for c=FNeg
        Regressors=[SAIF(MIdxs(c,1),:); squeeze(SHConvd(MIdxs(c,1),Keps1I(MIdxs(c,2)),:))'];
        CXs(:,c)=lsqnonneg(Regressors',DataToFit(c,:)');
    end
else
    for c=FNeg
        Regressors=[SAIF(MIdxs(c,1),:); squeeze(SHConvd(MIdxs(c,1),Keps1I(MIdxs(c,2)),:))'];
        CXs(:,c)=lsqnonneg(Regressors',DataToFit(c,:)',X(:,c)');
        %             X(:,c)=lsqnonneg(Regressors',DataToFit(c,:)',X(:,c)');
        %             X(:,c)=blocknnls(Regressors',DataToFit(c,:)');
        
        %             param = struct('A',Regressors','b',DataToFit(c,:)');
        %             [X(:,c),status] = asa_wrapper( X(:,c), zeros(2,1),
        %             Inf(2,1),'asa_quadratic_fcn','asa_quadratic_grad',
        %             'asa_quadratic_fcnGrad', struct('PrintParms',false), struct('PrintParms',false), param);
    end
end
Sims=zeros(nToFit,size(DataToFit,2));
if(nargout>2)
    for c=1:nToFit
        Regressors=[SAIF(MIdxs(c,1),:); squeeze(SHConvd(MIdxs(c,1),Keps1I(MIdxs(c,2)),:))'];
        Sims(c,:)=((Regressors')*CXs(:,c));
    end
end
warning(WarningS);

% if(ShowFig)
%     figure(1);clf;
%     nPlots=nToFit;
%     for i=1:nPlots
%         IStart=0;
%         gsubplot(nPlots,i);
%         plot(SampleTs,[DataToFit(i+IStart,:); Sims(i+IStart,:)]');
%         title([num2str(MIdxs(i+IStart,:))]);
%     end
% end
