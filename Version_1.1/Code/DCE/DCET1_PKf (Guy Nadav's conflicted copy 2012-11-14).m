PKFN=[WorkingP 'AfterPK' USStr '.mat'];
disp(PKFN);
% if(exist(PKFN,'file') && ~ForcePK)
%     D=dir(PKFN);
% %     if(D.datenum>DDateForPK)
% %         loadBut(PKFN,{CurVars.name});
% %         WorkingP=['/u/peptibase3-ext/libermg1/DCET1/' DCEInfos(CurSet).ShortStudyName filesep];
%         disp('PK loaded');
%         return;
% %     end
% end
% delete(PKFN);
% delete([WorkingP 'AfterCTC.mat']);

clear DCE4D pnk mixture DCE2D RDCE2D ExtendedMixture100 C1 C2 mixture100 GMixture optmixture100
disp(['DCE_PK, start ' datestr(now)]);

load(PrepareFN,'TimeBetweenDCEVols','BolusStart');
BolusStartSec=BolusStart*TimeBetweenDCEVols;
TimeBetweenDCEVols=TimeBetweenDCEVols*UnderSampling;
BolusStart=ceil(BolusStartSec/TimeBetweenDCEVols);
%% Simple noise measure
adCTC2D=abs(diff(CTC2D,[],2));
madCTC2D=median(adCTC2D,2);
rmadCTC2D=madCTC2D./max(CTC2D,[],2);
%%
% BolusStart=(BolusStart-1)*InterpolationFactor+1;
% 
% Base=3*InterpolationFactor+1+1;
% [a b]=max(SICTC2D(:,(BolusStart-3*InterpolationFactor-1):(BolusStart+3*InterpolationFactor+1)),[],2);
% N=histc(b,1:6*InterpolationFactor+2);
% [Tmp M]=max(N);
% BolusStart=BolusStart+M-Base;
% 
% Base=3*InterpolationFactor+1+1;
% [a b]=max(SICTC2D(:,(BolusStart-3*InterpolationFactor-1):(BolusStart+3*InterpolationFactor+1)),[],2);
% 
% Dif=b-Base;
% GoodRows=Dif>(-2*InterpolationFactor-1) & Dif<(2*InterpolationFactor+1) & Dif~=0;
% F=find(GoodRows);
% [U IA IB]=unique(Dif(GoodRows));
% for i=1:numel(U)
%     CDif=U(i);
%     CurIdxs=F(IB==i);
%     if(CDif<0)
%         SICTC2D(CurIdxs,(-CDif+1):end)=SICTC2D(CurIdxs,1:(end+CDif));
%         SICTC2D(CurIdxs,1:(-CDif))=0;
%     else
%         SICTC2D(CurIdxs,1:(end-CDif))=SICTC2D(CurIdxs,(CDif+1):end);
%         SICTC2D(CurIdxs,(end-CDif+1):end)=repmat(SICTC2D(CurIdxs,(end-CDif)),[1 CDif]);
%     end
% end
%%
% CTC2D=SICTC2D;
% nVols=size(SICTC2D,2);
% TimeBetweenDCEVols=TimeBetweenDCEVols/InterpolationFactor;
% BolusStartFromPrepare=BolusStart;
%%
nVols=size(CTC2D,2);
STimeBetweenDCEVols=TimeBetweenDCEVols;
STimeBetweenDCEVolsMin=TimeBetweenDCEVols/60;
SMin2Timestep=1/STimeBetweenDCEVolsMin;
nSVols=nVols;
SampleTs=((1:nSVols)-1)*STimeBetweenDCEVolsMin;
ShowFigures=false;
%% Basic definitions
Sdt=0.01;
SEndTime=10; % minutes
Sts=-SEndTime:Sdt:SEndTime; %(-nSVols*2:dVols:nSVols*2)./Min2Timestep;
options = struct('GradObj','off','Display','off');
InterestingTimePoints=Sts>-1;

Sts=SampleTs;

% A~mM.min, T~min; sig~min, alpha~mM, beta,s~1/min, tau~min
T1=0.17046;A1=0.809;sig1=0.0563;sig2=0.132;A2=0.330;T2=0.365;alpha=1.050;beta=0.1685;s=38.078;tau=0.483;
tauDelta=tau-T1;T2Delta=T2-T1;
% AIFParamsA=[T1 A1 sig1];
T1=1;
% AIF_Parker3=@(x) AIF_Parker(Sts,A1,x(3)*sig1,x(1)*T1,A2,x(3)*sig2,x(1)*T1+T2Delta,alpha,beta,s,x(1)*T1+tauDelta)*x(2)/1000;
% AIF_Parker6=@(x) AIF_Parker(Sts,A1,x(3)*sig1,x(1)*T1,A2*x(5),x(3)*sig2*x(6),x(1)*T1+T2Delta*x(4),alpha,beta,s,x(1)*T1+tauDelta)*x(2)/1000;
% AIF_Parker10=@(x) AIF_Parker( Sts,A1,x(3)*sig1,x(1)*T1,A2*x(5),x(3)*sig2*x(6),x(1)*T1+T2Delta*x(4),alpha*x(7),beta*x(8),s*x(9),x(1)*T1+tauDelta*x(10) )*x(2)/1000;
% OrigParkerAIF=AIF_Parker3(ones(1,3));

ParamAIFCoeff=[0.8,25,3;... % First bolus time
    0.25,4,1;... % First bolus magnitude
    0.25,10,2;... % First bolus width
    0.5,10,1;... % 2nd bolus time
    0.25,4,1;... % 2nd bolus magnitude
    0.25,4,1;... % 2nd bolus width
    0.5,2,1;... % alpha
    0.5,2,1;... % beta
    0.5,2,1;... % s
    0.5,5,1]; % tauDelta

ParamAIFCoeff(1,3)=BolusStart*TimeBetweenDCEVols/60;
ParamAIFCoeff(3,:)=[1 5 3];
ParamAIFCoeff(1,1)=0.1;
ParamAIFCoeff(1,2)=SampleTs(end)/2;
ParamAIFCoeff(5,1)=0;
ParamAIFCoeff(7,:)=[0.2 2 0.4];
ParamAIFCoeff(8,:)=[0 0.3 beta];
% ParamAIFCoeff(7,2)=20;
% ParamAIFCoeff(8,2)=20;
ParamAIFCoeff(10,1)=0;
ParamAIFCoeff(10,3)=0.5;
ParamAIFCoeff(2,2)=10;
ParamAIFCoeff(6,:)=[1 5 3];

nAIFParams=size(ParamAIFCoeff,1);
%% Options
% nSVox=nFClusters;
% SCWs=ones(1,nSVox)./nSVox;
% DesiredPseudoTimeBetweenDCEVols=1;
% nOptsKepT1=[21 22];
%%
Idx3D=NaN(size(Msk));
Idx3D(Msk2)=1:sumn(Msk2);
% for i=1:size(ChosenVoxels,1)
%     CVI(i)=Idx3D(ChosenVoxels(i,1),ChosenVoxels(i,2),ChosenVoxels(i,3));
% end
% CVI=CVI(isfinite(CVI));

% CVI=1:30;
% RCTCE100=-CTC2D(CVI,:);
%% Auto points
% MaxAroundBolus=max(CTC2D(:,(BolusStart-1-InterpolationFactor):(BolusStart+3*InterpolationFactor+1)),[],2);
MaxAroundBolus=max(CTC2D(:,(BolusStart-3):(BolusStart+3)),[],2);
MaxAtEnd=max(CTC2D(:,(end-5):end),[],2);
%%
S=sort(MaxAroundBolus);
N=size(CTC2D,1);
Vals=linspace(1,N,5);
Edges=S(floor(Vals));
[H,Bin] = histc(MaxAroundBolus,Edges);

S=sort(MaxAtEnd);
N=size(CTC2D,1);
Vals=linspace(1,N,5);
Edges=S(floor(Vals));
[H,Bin2] = histc(MaxAtEnd,Edges);

Clust=(Bin-1)*5+Bin2;
Reps=[];
Noises=[];
for i=1:max(Clust)
    CurI=find(Clust==i);
    Ns(i)=numel(CurI);
    if(numel(CurI)<100)
        continue;
    end
    [S Ord]=sort(rmadCTC2D(CurI));
    Reps=[Reps CurI(Ord(1:2)')'];
    Noises=[Noises S(1:2)'];
end
[S Ord]=sort(Noises);
Rep2D=CTC2D(Reps(Ord),:)./repmat(max(CTC2D(Reps(Ord),:),[],2),[1 size(CTC2D,2)]);
Good=find(sum(Rep2D(:,BolusStart+2:end)<-0.001,2)<3);
% Rep2D=Rep2D(Good,:);
CVI=Reps(Ord(Good));
RCTCE100=CTC2D(CVI,:);
%% Prepare first level
DesiredPseudoTimeBetweenDCEVols=STimeBetweenDCEVols;

NRCDCE=repMulti(RCTCE100,1./max(RCTCE100,[],2));
% NoisyCs=FindBizzare(mean(abs(diff(NRCDCE(:,(BolusStart+10):end),1,2)),2),1);
DiffStuff=diff(NRCDCE(:,(BolusStart+10):end),1,2);
if(isempty(DiffStuff))
    error('Diff empty in Noise calculation, DCET1_PKf !');
end
NoiseVals=mean(abs(DiffStuff),2);
[QQQ, optmixture] = GaussianMixture(NoiseVals, 3, 0,false);
[QQQ, BigNoiseGroupI]=max([optmixture.cluster.mu]);
NoisyCs=find(sum(optmixture.pnk(:,BigNoiseGroupI),2)>0.5);

LessOneCs=[]; %find(sum(RCTCE100(:,BolusStart+2:end)<1,2)>nVols*0.1);
BadCs=union(NoisyCs,LessOneCs);
BadCs=[];

PseudoSampleFac=ceil(STimeBetweenDCEVols./DesiredPseudoTimeBetweenDCEVols);
PseudoTimeBetweenDCEVols=STimeBetweenDCEVols/PseudoSampleFac;
PseudoTimeBetweenDCEVolsMin=PseudoTimeBetweenDCEVols/60;
PseudoMin2Timestep=1/PseudoTimeBetweenDCEVolsMin;
PseudoSampleTs=0:PseudoTimeBetweenDCEVolsMin:SampleTs(end);
nPSVols=numel(PseudoSampleTs);
AllTimePoints=1:nPSVols;

BolusStartFromPrepare=BolusStart;

SecLevelRes=2*4*2+1;
ParamMat=[0,2.2,5 11 SecLevelRes; 0.15,1,2, 11, SecLevelRes;0,0.3,1,17,SecLevelRes];
% BolusStart=floor(0.5*SMin2Timestep); % initial guess
BolusStart=BolusStartFromPrepare;


[SMaxR SMaxT]=max(RCTCE100,[],2);
% Rend=RCDCE100(:,end);
% REndRMaxR=Rend./SMaxR;
% figure;hist(SMaxR);
% BigMaxRVoxI=find(SMaxR>3);
[S Ord]=sort(SMaxR,'descend');
Ord=Ord(~ismember(Ord,BadCs))';
BigMaxRVoxISoon=Ord(SMaxT(Ord)<(nVols/2));
BigMaxRVoxI=union(Ord(1:min(5,numel(Ord))),BigMaxRVoxISoon(1:min(10,numel(BigMaxRVoxISoon))));
% BigMaxRVoxI=BigMaxRVoxISoon;
% BigMaxRVoxI=findTopIdxs(SMaxR,10);
% BigMaxRVoxI=getKrandomSamples(setdiff(BigMaxRVoxI,BadCs),min(numel(BigMaxRVoxI),5));
PreparedFullTimeAgT1=findPK_1LevelPrepareKnownT1(ParamMat,AllTimePoints,RCTCE100(BigMaxRVoxI,:),PseudoTimeBetweenDCEVols,[],1,PseudoSampleFac);