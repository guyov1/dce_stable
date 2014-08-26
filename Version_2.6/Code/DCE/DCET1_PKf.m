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

% Get the "derivative" of each pixel according to time (difference)
adCTC2D=abs(diff(CTC2D,[],2));
% Get the median of the derivative for each pixel in the entire time samples
madCTC2D=median(adCTC2D,2);
% For each median, calculate the percentage out of the maximal signal value during time of that pixel
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

%%Initial params
nVols=size(CTC2D,2);
STimeBetweenDCEVols=TimeBetweenDCEVols;
STimeBetweenDCEVolsMin=TimeBetweenDCEVols/60;
SMin2Timestep=1/STimeBetweenDCEVolsMin;
nSVols=nVols;
% Time sample in minutes for each temproal volume
SampleTs=((1:nSVols)-1)*STimeBetweenDCEVolsMin;
ShowFigures=false;

%% Basic definitions

% ASK GILAD - is Sdt supposed to be the time resolution? time intervals
%              between -10 to +10 minutes?
% ANSWER    - Yes.
Sdt=0.01;
SEndTime=10; % minutes
Sts=-SEndTime:Sdt:SEndTime; %(-nSVols*2:dVols:nSVols*2)./Min2Timestep;

options = struct('GradObj','off','Display','off');

% ASK GILAD - why did he choose ">-1" as the interesting time?
%             what is the meaning of negative time? The ability to shift
%             the curve?
% ANSWER    - He doesn’t really use it.
InterestingTimePoints=Sts>-1;
% Sts is now the time stamp in minutes for each temporal volume
Sts=SampleTs;

% A~mM.min, T~min; sig~min, alpha~mM, beta,s~1/min, tau~min
T1=0.17046;A1=0.809;sig1=0.0563;sig2=0.132;A2=0.330;T2=0.365;alpha=1.050;beta=0.1685;s=38.078;tau=0.483;
% ASK GILAD - What are those deltas? our T1,T2 against population average?
% ANSWER    - He is not sure if he used these specific variables exactly.
%             The T2-T1 is a better way to look on the second bolus than using absolute values.
%             The tau is for the sigmoid (which he doesn’t really use), "-T1" for a relative value and not absolute.
tauDelta=tau-T1;T2Delta=T2-T1;
% AIFParamsA=[T1 A1 sig1];

% ASK GILAD - Why did he assign T1=1?
% ANSWER    - T1 is the bolus start. He probably assumed it will be after a minute + he modifies it later.
T1=1;
% AIF_Parker3=@(x) AIF_Parker(Sts,A1,x(3)*sig1,x(1)*T1,A2,x(3)*sig2,x(1)*T1+T2Delta,alpha,beta,s,x(1)*T1+tauDelta)*x(2)/1000;
% AIF_Parker6=@(x) AIF_Parker(Sts,A1,x(3)*sig1,x(1)*T1,A2*x(5),x(3)*sig2*x(6),x(1)*T1+T2Delta*x(4),alpha,beta,s,x(1)*T1+tauDelta)*x(2)/1000;
% AIF_Parker10=@(x) AIF_Parker( Sts,A1,x(3)*sig1,x(1)*T1,A2*x(5),x(3)*sig2*x(6),x(1)*T1+T2Delta*x(4),alpha*x(7),beta*x(8),s*x(9),x(1)*T1+tauDelta*x(10) )*x(2)/1000;
% OrigParkerAIF=AIF_Parker3(ones(1,3));

% ASK GILAD - How did he determine those parameters? all of them, until
%             line 134...
% ANSWER    - Reasonable limits (lower bound, upper bound and start point).
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

% Updating the real bolus start time in minutes
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
% ASK GILAD - What "options" are these and what is the meaning of the
%             following?
% ANSWER    - He changed it in the newer version, so he doesn’t use it anymore.
DesiredPseudoTimeBetweenDCEVols=1;
nOptsKepT1=[21 22];


% NAN mask
Idx3D=NaN(size(Msk));
% Each pixel in the mask gets a linear index
Idx3D(Msk2)=1:sumn(Msk2);

% for i=1:size(ChosenVoxels,1)
%     CVI(i)=Idx3D(ChosenVoxels(i,1),ChosenVoxels(i,2),ChosenVoxels(i,3));
% end
% CVI=CVI(isfinite(CVI));

% CVI=1:30;
% RCTCE100=-CTC2D(CVI,:);

%% Get temporal interesting volumes around the bolus time and end time (it will characterize each C(t)).

% MaxAroundBolus=max(CTC2D(:,(BolusStart-1-InterpolationFactor):(BolusStart+3*InterpolationFactor+1)),[],2);
% ASK GILAD - What "Auto points" does he mean and why does he need the max 
%             of each pixel around the bolus temproal volumes and at the end? 
% ANSWER    - Each voxel has a c(t). The interesting characteristics of each voxel 
%             are the c(t) at bolus arrival and c(t) at the end of test (to know if we have enhancement).
%	          He used those values to cluster all the voxels so eventually he could work on some representatives.

MaxAroundBolus = max(CTC2D(:,(BolusStart-3):(BolusStart+3)),[],2);
MaxAtEnd       = max(CTC2D(:,(end-5):end),[],2);

%% Cluster each voxel according to the value around the bolus and end

% ASK GILAD - What was he trying to do in the following lines?
% ANSWER    - He tried to cluster the values around the bolus and at the end of test according to 5 bins.

% Sort the volumes around bolus time according to maximal value
S=sort(MaxAroundBolus);
% Number of masked pixels in each temproal volume
N=size(CTC2D,1);
% Line space 5 points between 1 and the number of masked pixels in each volume
Vals=linspace(1,N,5);
% Get the 5 values of the sorted maximal value pixels
Edges=S(floor(Vals));
% Histogram count for all the pixels around the bolus according to maximal
% and minimal value (1) divided to 5 bins
[H,Bin] = histc(MaxAroundBolus,Edges);

% Do the same sorting and binning for values close to the end of experiment
S        = sort(MaxAtEnd);
N        = size(CTC2D,1);
Vals     = linspace(1,N,5);
Edges    = S(floor(Vals));
[H,Bin2] = histc(MaxAtEnd,Edges);

% ASK GILAD - What is the meaning of the following clustering?
% ANSWER    - He gave a global clustering (both for bolus time and end time). 
%                       So if one point has index (2,3), it is actually 2*5 + 3.
Clust=(Bin-1)*5+Bin2;

Reps=[];
Noises=[];

%% Get representing voxels according to lowest noise

% For each cluster, get the least noisiest 2 voxels (and their index)
for i=1:max(Clust)
    
    % Indices of all voxels in cluster i
    CurI=find(Clust==i);
    
    % Number of voxels in that cluster
    Ns(i)=numel(CurI);
    
    % If there are not enough voxels in this cluster (<100), continue
    if(numel(CurI)<100)
        continue;
    end
    
    % Get a sorted list (S) of the noise medians of all the voxels from the
    % specific cluster. Get the original indexing (Ord) as well.
    [S Ord]=sort(rmadCTC2D(CurI));
    
    % Get the voxels indices with the loweset noise values (2 smallest)
    Reps=[Reps CurI(Ord(1:2)')'];
    
    % Get the noise values
    Noises=[Noises S(1:2)'];
end

% Sort the noises from all clusters (+ original indices)
[S Ord]=sort(Noises);


% Normalize each representing voxel in it's maximum value over time. 
% This way, their value will be between 0 -> 1
% Rep2D will hold the representing voxels over time (with normalized voxel
% values between 0 -> 1).
Rep2D=CTC2D(Reps(Ord),:)./repmat(max(CTC2D(Reps(Ord),:),[],2),[1 size(CTC2D,2)]);
% Good will hold all voxels that are not bad.
% Bad is defined as a voxel with negative C(t) in 3 or more time points.
Good=find(sum(Rep2D(:,BolusStart+2:end)<-0.001,2)<3);

% Rep2D=Rep2D(Good,:);
% Get the good voxels indices
CVI=Reps(Ord(Good));
% C(t) of all representing voxels
RCTCE100=CTC2D(CVI,:);

%% Prepare first level

% Time sample in seconds between each temproal volume
DesiredPseudoTimeBetweenDCEVols=STimeBetweenDCEVols;

% ASK GILAD - RCTCE100 is already normalized to values between 0->1. 
%                           Why is he normalizing again?
% ANSWER    - Seems like it is not normalized before. So here we normalize the C(t)s
NRCDCE=repMulti(RCTCE100, 1./max(RCTCE100,[],2));
% NoisyCs=FindBizzare(mean(abs(diff(NRCDCE(:,(BolusStart+10):end),1,2)),2),1);

% Get the derivative of C(t) 10 time points after the bolus to the end
DiffStuff=diff(NRCDCE(:,(BolusStart+10):end),1,2);

% if(isempty(DiffStuff))
%     error('Diff empty in Noise calculation, DCET1_PKf !');
% end

% Get the mean noise values after the bolus to the end
% NoiseVals=mean(abs(DiffStuff),2);

NoiseVals=EstimateNoise(RCTCE100);

% ASK GILAD - I didn't understand the following gaussians usage for the noise.
% ANSWER -    He used it for clustering the noise (to work on a smaller group)

% Create a gaussian mixture for the noise
[QQQ, optmixture] = GaussianMixture(NoiseVals, 3, 0,false);

% Get the biggest noise group (if there is more than one)
% QQQ will hold the noise value and BigNoiseGroupI the index
[QQQ, BigNoiseGroupI]=max([optmixture.cluster.mu]);

% Find the noisy C(t)s
% Should uncomment when getting the latest GaussianMixture from Gilad
% Currentyly, the "pnk" does not exist
%NoisyCs=find(sum(optmixture.pnk(:,BigNoiseGroupI),2)>0.5);

% ASK GILAD - What did he try to do here? Seems like he did not pick any bad C(t)s...
% ANSWER - He didnt use it at the end.
LessOneCs=[]; %find(sum(RCTCE100(:,BolusStart+2:end)<1,2)>nVols*0.1);

% Should uncomment later
%BadCs=union(NoisyCs,LessOneCs);
BadCs=[];

% ASK GILAD - What is the meaning of the following var if he defined both
%                           DesiredPseudoTimeBetweenDCEVols and STimeBetweenDCEVols to get the same value.
% ANSWER - He didnt use it at the end.
PseudoSampleFac=ceil(STimeBetweenDCEVols./DesiredPseudoTimeBetweenDCEVols);
PseudoTimeBetweenDCEVols=STimeBetweenDCEVols/PseudoSampleFac;
PseudoTimeBetweenDCEVolsMin=PseudoTimeBetweenDCEVols/60;
% Number of temporal samples in 1 minute
PseudoMin2Timestep=1/PseudoTimeBetweenDCEVolsMin;
% All time samples of C(t)s - currently the same as the test itself
PseudoSampleTs=0:PseudoTimeBetweenDCEVolsMin:SampleTs(end);
% Number of time samples
nPSVols=numel(PseudoSampleTs);
% Index of each time point
AllTimePoints=1:nPSVols;
% The temporal volume of bolus start than we prepared before
BolusStartFromPrepare=BolusStart;
% ASK GILAD - What is the meaning of the following parameters?
% ANSWER - Gilad says this is old garbage. He chaned it after using Murase.
SecLevelRes=2*4*2+1;
ParamMat=[0,2.2,5 11 SecLevelRes; 0.15,1,2, 11, SecLevelRes;0,0.3,1,17,SecLevelRes];
% BolusStart=floor(0.5*SMin2Timestep); % initial guess
% ASK GILAD - Is this redundant? He just calculated the opposite a few lines above...
% ANSWER - Yes.
BolusStart=BolusStartFromPrepare;

% Get the max of each C(t) of the representing voxels (and the index)
[SMaxR SMaxT]=max(RCTCE100,[],2);
% Rend=RCDCE100(:,end);
% REndRMaxR=Rend./SMaxR;
% figure;hist(SMaxR);
% BigMaxRVoxI=find(SMaxR>3);

% Sort the maximal values of all C(t)s and keep their original index
[S Ord]=sort(SMaxR,'descend');
% Screen the bad C(t)s
Ord=Ord(~ismember(Ord,BadCs))';

% ASK GILAD - What is the meaning of the next 2 lines?
% ANSWER - He tried to take the voxels with high values (artery)
BigMaxRVoxISoon=Ord(SMaxT(Ord)<(nVols/2));
BigMaxRVoxI=union(Ord(1:min(5,numel(Ord))),BigMaxRVoxISoon(1:min(10,numel(BigMaxRVoxISoon))));

% BigMaxRVoxI=BigMaxRVoxISoon;
% BigMaxRVoxI=findTopIdxs(SMaxR,10);
% BigMaxRVoxI=getKrandomSamples(setdiff(BigMaxRVoxI,BadCs),min(numel(BigMaxRVoxI),5));

%% Gilad says this is old and only calculates a few variables that will be used later.
PreparedFullTimeAgT1=findPK_1LevelPrepareKnownT1(ParamMat,AllTimePoints,RCTCE100(BigMaxRVoxI,:),PseudoTimeBetweenDCEVols,[],1,PseudoSampleFac);