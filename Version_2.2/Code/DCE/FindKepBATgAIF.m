function [MIdxs CXs Sims] = FindKepBATgAIF(DataToFit,SAIF,SHConvd,Keps1I)

% BAT Kep Vp Ktrans Calculation

%% Initiate variables 

% Number of representing voxels to fit
nToFit=size(DataToFit,1);
% Number of possible keps
nKeps1=size(SHConvd,2);
% Number of possible shift in time (delta T's)
nTDif=size(SAIF,1);
% RMSs will hold the root squared error for each time shift and kep
RMSs=zeros([nTDif nKeps1 nToFit]);
% Xs will hold te paramters for each time shift(i) and kep(j)
Xs=zeros([nTDif nKeps1 2 nToFit]);

WarningS=warning('off','MATLAB:rankDeficientMatrix');

%% Calculate Vp and Ktrans for each time shift and kep

% Go over all possible delta T's
for i=1:nTDif
    
    %     disp([i nTDif]);
    % Duplicate the shifted Parker's AIFs
    Regressors=[SAIF(i,:); SAIF(i,:)];
    
    % For each possible kep
    for j=1:nKeps1
        % Put the exp(-kep*t) CONV AIF(t)  with the specific delta t (i)
        % and kep (j) in the second Regressors line
        Regressors(2,:)=SHConvd(i,Keps1I(j),:);

        % ASK GILAD - Make sure the following explanation I wrote is what he meant.
        % ASNWER - Yes.
        % I think Gilad took the follwing equation from Fluckiger:
        % Ct(t) = Ktrans*Cp(t) CONV exp(-kep*t) + VpCp(t) 
        % Ct(t) is supposed to be the representing voxels. Cp(t) the AIF.
        %
        % DataToFit=X' * Regressors
        %  
        %  [DataToFit ]   =  [ Vp     Ktrans    ]      *     [            AIF(t)                                            ]
        %  [        ...              ]   =  [  ..             ..            ]             [             AIF(t) CONV  exp(-kep*t)      ]        
        %
        %The following gives us:
        %    X' =   [ Vp                ... ]
        %              [ Ktrans      ... ]
        X=Regressors'\DataToFit';
        
        % Xs will hold te paramters for each time shift(i) and kep(j)
        Xs(i,j,:,:)=X;
        % Calculated the diff and RMS of the result compating to the real DataToFit
        Difs=((Regressors')*X)'-DataToFit;
        RMSs(i,j,:)=sqrt(mean(Difs.^2,2));
    end
    
end

%% Pick vp,ktrans with lowest RMS for each representing voxel and make sure
%% it is non-negative

% CXs will hold for each representing voxel the best vp and ktrans (minimal RMS)
CXs=zeros([2 nToFit ]);
% MIdxs will hold the indices of the best vp,ktrans out of all possible delta t's and kep's
MIdxs=zeros(nToFit,2);

% For each representing voxel
for c=1:nToFit
    % gmin returns the minimum RMS (and index) over all delta T's and kep's
    % ( we do it for each representing voxel number "c")
    [Tmp MIdxs(c,:)]=gmin(squeeze(RMSs(:,:,c)));
    % For the current representing voxel, take the vp and ktrans with minimal RMS
    CXs(:,c)=squeeze(Xs(MIdxs(c,1),MIdxs(c,2),:,c));
end

% Check if any of the columns (of the representing voxels) has a negative value
FNeg=find(any(CXs<0,1));
v=ver('MATLAB');
% For each negative parameter we got, calculate least squares that minimize
% the parameters under the constraint they have to be positive
if(str2num(v.Release(5:6))>10)
    for c=FNeg
		% Get the relevant Regressors for the calculation
        Regressors=[SAIF(MIdxs(c,1),:); squeeze(SHConvd(MIdxs(c,1),Keps1I(MIdxs(c,2)),:))'];
        CXs(:,c)=lsqnonneg(Regressors',DataToFit(c,:)');
    end
else
    for c=FNeg
		% Get the relevant Regressors for the calculation
        Regressors=[SAIF(MIdxs(c,1),:); squeeze(SHConvd(MIdxs(c,1),Keps1I(MIdxs(c,2)),:))'];
		% ASK GILAD - 1. Is least squares equivalent to the method above where
	    %                                he used X=Regressors'\DataToFit' ?
	    %                          2. I dont understand the starting point it gave. Didnt he mean to give it the negative
	    %                              parameters which are the optimal? Currently he gives "X" which holds the last time shift
	    %                              and kep calculation from before (and is not related to the negative value).
	    % ANSWER - 1. Yea. Dividing is like least squares just that it gives negative values.
	    %                    2. He ment to give the relevant negative kep. Should fix it if it is not behaving this way...
	    %
	    % Get new parametrs instead of the negative ones using least squares non-negative
        CXs(:,c)=lsqnonneg(Regressors',DataToFit(c,:)',X(:,c)');
        %             X(:,c)=lsqnonneg(Regressors',DataToFit(c,:)',X(:,c)');
        %             X(:,c)=blocknnls(Regressors',DataToFit(c,:)');
        
        %             param = struct('A',Regressors','b',DataToFit(c,:)');
        %             [X(:,c),status] = asa_wrapper( X(:,c), zeros(2,1),
        %             Inf(2,1),'asa_quadratic_fcn','asa_quadratic_grad',
        %             'asa_quadratic_fcnGrad', struct('PrintParms',false), struct('PrintParms',false), param);
    end
end

%% Calculare the new AIF for each representing voxel out of the calculted vp,ktrans

%  Sims will hold the fitted AIF using the calculated parameters for each representing voxel 
Sims=zeros(nToFit,size(DataToFit,2));

% Sims is the third output argument
if(nargout>2)
    
    % For each representing voxel
    for c=1:nToFit
        
        % Sims will hold the fitted AIF using the calculated parameters for
        % each representing voxel
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
