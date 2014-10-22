function WholeVolCompute_ManualArtNoBAT(WorkingP)
% WorkingP=['/data/home/gilad/Desktop/DCE/gilad/Database/DCEOut/' D(pp).name filesep];
PKMFN=[WorkingP 'ManualArtNoBATPK2.mat'];
if(exist(PKMFN,'file'))
    disp([D(pp).name ' already exist']);
    return;
end
if(~exist([WorkingP 'AfterCTC.mat'],'file'))
    error('No AfterCTC');
    return;
end
a=load([WorkingP 'AfterCTC.mat']);

nVols=size(a.CTC2D,2);
STimeBetweenDCEVols=a.TimeBetweenDCEVols;
STimeBetweenDCEVolsMin=a.TimeBetweenDCEVols/60;
SMin2Timestep=1/STimeBetweenDCEVolsMin;
nSVols=nVols;
% if(nSVols>70 || nSVols<30)
%     continue;
% end
SampleTs=((1:nSVols)-1)*STimeBetweenDCEVolsMin;

PKMFNAIF=[WorkingP 'ManualArtNoBATPK_AIF.mat'];
load(PKMFNAIF);
% MeanArtCTC=MeanArtCTC2;

N=size(a.CTC2DBigGood,1);
PKs=NaN(N,24);

CSAIF=cumtrapz(SampleTs,MeanArtCTC);

MskX=a.DBrainMask;
MskX(MskX)=a.MskCTCGood;
%

NAtATime=5000;
disp(['There are ' num2str(N) ' voxels to compute ' PKMFN]);
before=now;


% for i=1:NAtATime:N
%     tic
%     CurIs=i:min(N,i+NAtATime-1);
%     PKs(CurIs,:) = FindPKBATgAIFMuraseF4Models_TProb(a.CTC2DBigGood(CurIs,:),MeanArtCTC,SampleTs,CSAIF);
%     t=toc;
%     TimeFromStart=now-before;
%     WillFinishAt=before+TimeFromStart*N/CurIs(end);
%     disp(['Calculating ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
% end

%% --------------------------- ParFor Try! --------------------------------
num_iter = floor(N/NAtATime);
PKs_check = zeros(num_iter,NAtATime,size(PKs,2));
CTC2DBigGood_Sliced = zeros(num_iter,NAtATime,size(a.CTC2DBigGood,2));
% Copy CTC2DBigGood into a sliced variable
for k=1:num_iter
    index = 1 + ( (k-1)*NAtATime );
    CurIs=index:index+NAtATime-1;
    CTC2DBigGood_Sliced(k,:,:) = a.CTC2DBigGood(CurIs,:);
end

parfor j=1:num_iter
    tic
    index = 1 + ( (j-1)*NAtATime );
    CurIs=index:index+NAtATime-1;
    
    PKs_check(j,:,:) = FindPKBATgAIFMuraseF4Models_TProb(squeeze(CTC2DBigGood_Sliced(j,:,:)),MeanArtCTC,SampleTs,CSAIF);
    %PKs(CurIs,:) = FindPKBATgAIFMuraseF4Models_TProb(CTC2DBigGood(CurIs,:),SAIF,SampleTs,CSAIF);
    
    t=toc;
    TimeFromStart=now-before;
    WillFinishAt=before+TimeFromStart*N/CurIs(end);
    disp(['Calculating manual BAT ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
end

% Put the resulting paramters in the original variable
for k=1:num_iter
    index = 1 + ( (k-1)*NAtATime );
    CurIs=index:index+NAtATime-1;
    PKs(CurIs,:) = PKs_check(k,:,:);
end

% Handle last iteration
last_index = N - mod(N,NAtATime) + 1;
tic;
CurIs=last_index:N;
PKs(CurIs,:) = FindPKBATgAIFMuraseF4Models_TProb(a.CTC2DBigGood(CurIs,:),SAIF,SampleTs,CSAIF);
t=toc;
TimeFromStart=now-before;
WillFinishAt=before+TimeFromStart*N/CurIs(end);
disp(['Calculating manual BAT ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);

%% -----------------------------------------------------------------------------------------------


disp(['Finished ' WorkingP]);
save(PKMFN,'PKs','MskX','MeanArtCTC');
disp('Computation finshed, making niis...');
%%
TTls={'BATfinal' 'VpFinal' 'KtransFinal' 'Kepfinal' 'VeFinal' 'RSSFinal' 'RSS0' 'RSS1' 'RSS2' 'RSS3' 'F1v0' 'F2v1' 'F3v2' 'BAT1' 'Vp1' 'BAT2' 'Vp2' 'Ktrans2' 'BAT3' 'Vp3' 'Ktrans3' 'Kep3' 'Ve3' 'WhichModel'};

RelaxP=[WorkingP 'ManualArtNoBAT2' filesep];
mkdir(RelaxP);
MeanFN=[WorkingP 'DCEMean.nii'];
disp(RelaxP);

N=sumn(MskX);
for i=1:numel(TTls)
    Tmp3D=MskX*0;
    Tmp3D(MskX)=PKs(1:N,i);
    Raw2Nii(Tmp3D,[RelaxP TTls{i} '.nii'],'float32', MeanFN);
end
disp('Finished Niftis');