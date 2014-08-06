function FindAIFByManualArtNoBAT(WorkingP,DoFull3D)

% WorkingP=['/data/home/gilad/Desktop/DCE/gilad/Database/DCEOut/' D(pp).name filesep];
PKMFN=[WorkingP 'ManualArtNoBATPK_AIF.mat'];
if(exist(PKMFN,'file') && ~DoFull3D)
%     disp([PKMFN ' already exist']);
%     return;
end

After_CTC_file = dir([WorkingP 'AfterCTC*.mat']);
After_CTC_file = After_CTC_file.name;

%if(~exist([WorkingP 'AfterCTC.mat'],'file'))
if(~exist([WorkingP After_CTC_file],'file'))
    error('No AfterCTC');
    return;
end

%a=load([WorkingP 'AfterCTC.mat']);
a=load([WorkingP After_CTC_file]);


nVols=size(a.CTC2D,2);
STimeBetweenDCEVols=a.TimeBetweenDCEVols;
STimeBetweenDCEVolsMin=a.TimeBetweenDCEVols/60;
SMin2Timestep=1/STimeBetweenDCEVolsMin;
nSVols=nVols;
% if(nSVols>70 || nSVols<30)
%     continue;
% end
SampleTs=((1:nSVols)-1)*STimeBetweenDCEVolsMin;

ManualArt3D=loadniidata([WorkingP 'manualArt.nii']);
F1=find(ManualArt3D);
F2=find(a.Msk2);
CVI=find(ismember(F2,F1));
ArtCTCs=a.CTC2D(CVI,:);
MeanArtCTC=mean(ArtCTCs,1);
MeanArtCTC=MeanArtCTC*max(ArtCTCs(:))./max(MeanArtCTC);
MeanArtCTC2=mean(NormalizeByRows(ArtCTCs),1);
MeanArtCTC2=MeanArtCTC2*max(ArtCTCs(:))./max(MeanArtCTC2);

figure(100);clf;subplot(2,2,1);plot(SampleTs,ArtCTCs');subplot(2,2,2);plot(SampleTs,NormalizeByRows(ArtCTCs)');
subplot(2,2,3);plot(SampleTs,ArtCTCs');hold on;plot(SampleTs,MeanArtCTC,'k','LineWidth',2);plot(SampleTs,MeanArtCTC2,'m','LineWidth',2);
subplot(2,2,4);plot(SampleTs,NormalizeByRows(ArtCTCs)'*max(ArtCTCs(:)));hold on;plot(SampleTs,MeanArtCTC,'k','LineWidth',2);plot(SampleTs,MeanArtCTC2,'m','LineWidth',2);
gprint(100,[WorkingP 'manualArt.png']);
close(100);
MeanArtCTC=MeanArtCTC2;
save(PKMFN,'MeanArtCTC');
if(~DoFull3D)
    return;
end

N=size(a.CTC2DBigGood,1);
PKs=NaN(N,28);

CSAIF=cumtrapz(SampleTs,MeanArtCTC);

MskX=a.DBrainMask;
MskX(MskX)=a.MskCTCGood;
%

NAtATime=5000;
disp(['There are ' num2str(N) ' voxels to compute ' PKMFN]);
before=now;
for i=1:NAtATime:N
    tic
    CurIs=i:min(N,i+NAtATime-1);
    PKs(CurIs,:) = FindPKBATgAIFMuraseF4Models_TProb(a.CTC2DBigGood(CurIs,:),MeanArtCTC,SampleTs,CSAIF);
    t=toc;
    TimeFromStart=now-before;
    WillFinishAt=before+TimeFromStart*N/CurIs(end);
    disp(['Calculating Manual No BAT ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
end
disp(['Finished ' WorkingP]);
save(PKMFN,'PKs','MskX','MeanArtCTC');

TTls={'BATfinal' 'VpFinal' 'KtransFinal' 'Kepfinal' 'VeFinal' 'RSSFinal' 'RSS0' 'RSS1' 'RSS2' 'RSS3' 'F1v0' 'F2v1' 'F3v2' 'BAT1' 'Vp1' 'BAT2' 'Vp2' 'Ktrans2' 'BAT3' 'Vp3' 'Ktrans3' 'Kep3' 'Ve3' 'WhichModel'};
RelaxP=[WorkingP 'ManualArtNoBAT3' filesep];
mkdir(RelaxP);
MeanFN=[WorkingP 'DCEMean.nii'];
disp(RelaxP);

a=load(PKMFN);
N=sumn(MskX);
for i=1:numel(TTls)
    Tmp3D=MskX*0;
    Tmp3D(MskX)=PKs(1:N,i);
    Raw2Nii(Tmp3D,[RelaxP TTls{i} '.nii'],'float32', MeanFN);
end
disp('Finished Manual No BAT Niftis');