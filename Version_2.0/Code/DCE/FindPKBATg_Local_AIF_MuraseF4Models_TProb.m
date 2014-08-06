function [PKs, Sims, Cost] = FindPKBATg_Local_AIF_MuraseF4Models_TProb(DataToFit,SAIF,SampleTs,CumSumAIFdTSamp,WPerCTC,MTT_map)
                           
%% Find PK parameters using diffrence of C, as in Murase 2004

% Disable rank Deficient warning
WarningS=warning('off','MATLAB:rankDeficientMatrix');

% dT=SampleTs(2)-SampleTs(1);

nToFit     = size(DataToFit,1);
nVols      = size(DataToFit,2);
nTDif      = size(SAIF,1);
BADInModel = nTDif>0;
alpha      = [0.01 0.01 0.000001];
df1        = 1+BADInModel;
df2        = 2+BADInModel;
df3        = 3+BADInModel;
Denom1     = nVols-df1;
Denom2     = nVols-df2;
Denom3     = nVols-df3;

%    1        2         3        4          5     6        7    8   9    10
% BATfinal VpFinal KtransFinal Kepfinal VeFinal RSSFinal RSS0 RSS1 RSS2 RSS3
%  11   12   13   14  15   16   17   18     19  20      21   22  23   24
% F1v0 F2v1 F3v2 BAT1 Vp1 BAT2 Vp2 Ktrans2 BAT3 Vp3 Ktrans3 Kep3 Ve3 WhichModel
%   25     26      27        28
% TProb1 TProb2  TProb3  TProbFinal

PKs  = NaN(nToFit,28);
Sims = zeros(size(DataToFit));

Vp1     = NaN(1,nTDif);
Vp2     = NaN(1,nTDif);
Vp3     = NaN(1,nTDif);
RSS1    = NaN(1,nTDif);
RSS2    = NaN(1,nTDif);
RSS3    = NaN(1,nTDif);
Ktrans2 = NaN(1,nTDif);
Ktrans3 = NaN(1,nTDif);
Kep3    = NaN(1,nTDif);

% mCumSumDatadT=cumtrapz(DataToFit,2).*dT;
mCumSumDatadT = -cumtrapz(SampleTs,DataToFit')';
% mCumSumDatadT=-[zeros(nToFit,1) cumsum(DataToFit(:,1:end-1),2)]; % *dT;
SimsModels    = zeros(4,nVols);

% Go over every voxel
for c=1:nToFit
    
    CurData          = DataToFit(c,:)';
    CurmCumSumDatadT = mCumSumDatadT(c,:)';
    
    % Model 0
    PKs(c,7) = sum(DataToFit(c,:).^2); % RSS0
    
    for t=1:nTDif
        CurMat=[CumSumAIFdTSamp(t,:)' CurmCumSumDatadT SAIF(t,:)'];
        % Model 1
        Vp1(t)  = max(0,CurMat(:,3)\CurData);
        RSS1(t) = sum((CurData-CurMat(:,3)*Vp1(t)).^2);
        
        % Model 2
        B=max(0,CurMat(:,[1 3])\CurData);
        Ktrans2(t)=B(1);
        Vp2(t)=B(2);
        RSS2(t)=sum((CurData-CurMat(:,3)*Vp2(t)-CurMat(:,1)*Ktrans2(t)).^2);
        
        % Model 3
        B          = CurMat\CurData;
        Vp3(t)     = max(0,B(3));
        Kep3(t)    = max(0,B(2));
        Ktrans3(t) = max(0,B(1)-Kep3(t)*Vp3(t));
        %         figure;plot(SampleTs,A(:,1),'k',SampleTs,A(:,2),'r')
        %         figure;plot(SampleTs,A(:,1),'k',SampleTs,A(:,2),'b',SampleTs,DataToFit(c,:),'r',SampleTs,SAIF(t,:),'g')
        RSS3(t)    = sum((CurData-CurMat(:,3)*Vp3(t)-CurMat(:,1)*(Ktrans3(t)+Vp3(t)*Kep3(t))-CurMat(:,2)*Kep3(t)).^2);
        %         figure;plot(SampleTs,DataToFit(c,:),'k',SampleTs,CForCluster(t,:),'b');
        %         figure;plot(SampleTs,[DataToFit(c,:); Ax(:,1)'*B(1); Ax(:,2)'*B(2); Ax(:,3)'*B(3);CForCluster(t,:)])
        % figure;plot(SampleTs,CurData,'k',SampleTs,CurMat(:,3)*Vp3(t),'r',SampleTs
        % ,CurMat(:,1)*(Ktrans3(t)+Vp3(t)*Kep3(t)),'g',SampleTs,CurMat(:,2)*Kep3(t)
        % ,'b',SampleTs,CurMat(:,3)*Vp3(t)+CurMat(:,1)*(Ktrans3(t)+Vp3(t)*Kep3(t))+CurMat(:,2)*Kep3(t),'m',SampleTs,CurMat(:,3)*Vp2(t)+CurMat(:,1)*Ktrans2(t),'c')
    end
    
    [PKs(c,8) PKs(c,14)] = min(RSS1); % RSS1 BAT1
    [PKs(c,9) PKs(c,16)] = min(RSS2); % RSS2 BAT2
    [PKs(c,10) PKs(c,19)]= min(RSS3); % RSS3 BAT3
    
    PKs(c,11)= fcdf(((PKs(c,7)-PKs(c,8))/df1)/(PKs(c,8)/Denom1),df1,Denom1); % F1v0
    PKs(c,12)= fcdf(((PKs(c,8)-PKs(c,9)))/(PKs(c,9)/Denom2),1,Denom2); % F2v1
    PKs(c,13)= fcdf(((PKs(c,9)-PKs(c,10)))/(PKs(c,10)/Denom3),1,Denom3); % F3v2
    
    WhichModel = find(PKs(c,11:13)<(1-alpha),1,'first'); % WhichModel
    if(isempty(WhichModel))
        WhichModel=4;
    end
    PKs(c,24) = WhichModel;
    
    PKs(c,15) = Vp1(PKs(c,14)); %Vp1
    PKs(c,17) = Vp2(PKs(c,16)); %Vp2
    PKs(c,20) = Vp3(PKs(c,19)); %Vp3
    PKs(c,18) = Ktrans2(PKs(c,16));% Ktrans2
    PKs(c,21) = Ktrans3(PKs(c,19));% Ktrans3
    PKs(c,22) = Kep3(PKs(c,19));% Kep3
    PKs(c,23) = PKs(c,21)/PKs(c,22);% Ve3
    
    SimsModels(2,:) = SAIF(PKs(c,14),:)*PKs(c,15);
    SimsModels(3,:) = SAIF(PKs(c,16),:)*PKs(c,17)+CumSumAIFdTSamp(PKs(c,16),:)*PKs(c,18);
    SimsModels(4,:) = SAIF(PKs(c,19),:)*PKs(c,20)+CumSumAIFdTSamp(PKs(c,19),:)*(PKs(c,21)+PKs(c,20)*PKs(c,22))+CurmCumSumDatadT'*PKs(c,22);
    
    BATs     = [NaN PKs(c,[14 16 19])];
    PKs(c,1) = BATs(WhichModel); % BATFinal
    RSSs     = PKs(c,7:10);
    PKs(c,6) = RSSs(WhichModel); % RSSFinal
    VPs      = [NaN PKs(c,[15 17 20])];
    PKs(c,2) = VPs(WhichModel); % VpFinal
    KTranss  = [NaN NaN PKs(c,[18 21])];
    PKs(c,3) = KTranss(WhichModel); % KtransFinal
    
    if(WhichModel==4)
        PKs(c,4)=PKs(c,22); % Kepfinal
        PKs(c,5)=PKs(c,3)./PKs(c,4); % VeFinal
    end
    
    %% TProb
    Others=setdiff(1:nTDif,(PKs(c,14)-2:PKs(c,14)+2));
    if(~isempty(Others))
        PKs(c,25)=PKs(c,8)./min(RSS1(Others)); % TProb1
    end
    
    Others=setdiff(1:nTDif,(PKs(c,16)-2:PKs(c,16)+2));
    if(~isempty(Others))
        PKs(c,26)=PKs(c,9)./min(RSS2(Others)); % TProb2
    end
    
    Others=setdiff(1:nTDif,(PKs(c,19)-2:PKs(c,19)+2));
    if(~isempty(Others))
        PKs(c,27)=PKs(c,10)./min(RSS3(Others)); % TProb3
    end
    TProbs=[NaN PKs(c,25:27)];
    PKs(c,28)=TProbs(WhichModel); % TProbFinal
    
    Sims(c,:)=SimsModels(WhichModel,:);
    %     figure;plot([CurData'; SimsModels]')
end

% Enable back rank Deficient warning
warning(WarningS);

% Calculate the RMS fit of the data to the simulated curve (Sims)
if(nargout>=3)
    Cost=max(sqrt(mean(((Sims-DataToFit).^2),2)).*WPerCTC);
end

% DEBUG STUFF
if(false)
    
    figure(111);
    clf;
    nPlots=nToFit;
    for i=1:nPlots
        IStart=0;
        gsubplot(nPlots,i);
        plot(SampleTs,DataToFit(i+IStart,:),'k','LineWidth',2);hold on;plot(SampleTs, Sims(i+IStart,:),'b');
        title([num2str(PKs(i+IStart,1))]);
    end
end