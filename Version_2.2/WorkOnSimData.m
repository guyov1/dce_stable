% StudyName='RoAs_20080122';
% StudyName='WhSa_20070813';
BasePath=[BaseBaseP StudyName filesep];
WorkingP=BasePath;
FinalT1=load([BasePath 'AfterCTC.mat'],'FinalT1');
FinalT1=FinalT1.FinalT1;
MskC=load([BasePath 'MskC.mat']);
MskC=MskC.MskC;
Sims=load([BasePath 'SimSigNoised.mat']);
Sims=Sims.Sims;
disp('Loaded simulated signal');
%%
% Sig to T1
SDCEFA=25;
SDCETR=5.6;
FA1D=FinalT1(MskC)*0+SDCEFA;
SM0=10000;
M01D=FinalT1(MskC)*0+SM0;

DCEVolToT1=@(FAx,M0x,Vol) SDCETR./(-log(((Vol./M0x)-sind(FAx))./((Vol./M0x).*cosd(FAx)-sind(FAx))));
T12D=zeros(size(Sims));
for i=1:size(Sims,2)
    disp([1 i]);
    T12D(:,i)=DCEVolToT1(FA1D,M01D,Sims(:,i));
end
% T1 to CTC
CTC2D=1./T12D-1./repmat(FinalT1(MskC),1,size(Sims,2));
clear T12D Sims
disp('Got CTC2D');
% q=10000;gfig;plot(CTC2D(q,:),'*');
%% Sample in temporal domain
HdtBase=0.5/60;
ReqDt=6/60;
SampleRatio=ReqDt/HdtBase;
CTC2D=CTC2D(:,1:SampleRatio:end);
save([BasePath 'TSampledSimCTCs.mat'],'CTC2D');
disp('Sampled in temporal domain');
%% Find arteries
%% Filters
Msk3D=MskC;
FBrainMask=bwfillHoles3Dby2D(MskC);
for i=1:size(FBrainMask,3)
    MskE(:,:,i)=imerode(squeeze(FBrainMask(:,:,i)),strel('disk',16));
end
NumVols=size(CTC2D,2);
TimeBetweenDCEVolsFinal=ReqDt*60;
TimeBetweenDCEVolsMin=ReqDt;
MskNotInEdge=MskE & Msk3D;

MedTC=median(CTC2D,1);
TwoMinTimePoint=floor(2/TimeBetweenDCEVolsMin);
Ps=zeros(1,numel(MedTC))+2;
% We use the t-test to get the biggest probability that the distribution of the sample
% is diffrent than the rest of the test ( -> smallest Ps value)
for i=3:min(TwoMinTimePoint,numel(MedTC)-2) %Take the minimum out of 2 minutes frame to 2 frames before the end
    [h Ps(i)]=ttest2(MedTC(1:i),MedTC((i+1):end),[],[],'unequal');
end
mLPs=-log(Ps);
% figure;plot(1:numel(MedTC),MedTC,'b',1:numel(MedTC),mLPs.*(max(MedTC)-min(MedTC))./(max(mLPs)-min(mLPs))+min(MedTC),'r')
[Tmp, BolusStart]=max(mLPs);


[DataOutA, ~, OtherMasks, Msks]=FilterCTCToFindArt(CTC2D(:,1:NumVols),Msk3D,NumVols,TimeBetweenDCEVolsFinal,BolusStart,struct('TimeDelayToMaskVeins',-0.5));
for i=1:numel(Msks)
    Msks{i}=Msks{i} & MskNotInEdge;
    MsksN(i)= sumn(Msks{i});
end
AllMsks=cat(4,Msks{:});
SumMsk=sum(AllMsks,4);
NinNsks=histc(SumMsk(:),0:numel(Msks));
TrgN=300;
if(NinNsks(end)>TrgN)
    MskF=SumMsk==numel(Msks);
else
    MskAI=numel(Msks)-find(cumsum(flipud(NinNsks))>TrgN,1)+1;
    MskA=SumMsk>MskAI;
    MskB=(SumMsk==MskAI);
    TmpF=find(MskB);
    Needed=TrgN-sumn(MskA);
    RIdxs=randsample(numel(TmpF),Needed);
    Tmp=TmpF*0>1;
    Tmp(RIdxs)=true;
    MskB(MskB)=Tmp;
    MskF=MskA | MskB;
end
% MaskWithEnoughI=find(MsksN>10,1,'last');
% MskF=Msks{MaskWithEnoughI} & MskNotInEdge;
% MskF=OtherMasks & MskNotInEdge;
InIdxs=getIndicesOfMskInsideMsk(MskF,Msk3D);
DataOut=CTC2D(InIdxs,:);
% Raw2Nii(MskF,[WorkingP 'ArtMsk.nii'],'float32',MeanFN);
% gfig;plot(NormalizeByRows(DataOut(200:300:end,:))')
disp('Calculated filters');
%%  Get the representing voxels data and noise
Defaults.SubSampling=1;
Defaults.nVolsToRemoveFromEnd=0;
Defaults.SubSecondResolution=2;
Defaults.MinFirstBolusStd=2;
Defaults.EM_Num_Of_Iterations=0;
Defaults.FMS_TolFun=1e-11;
Defaults.FMS_MaxFunEvals=10000;
Defaults.FMS_MaxIter=10000;
Defaults.MaxTDif_ForAIFSearch=3;
Defaults.MaxTDif_ForWholeVOI=6;
Defaults.Rep_MaxAroundBolus=4;
Defaults.Rep_RatioToEnd=4;
Defaults.Rep_nPerSet=2;
Defaults.Mask_thresh=0.5;
Defaults.Run_On_All=0;
Defaults.TimeDelayToMaskVeins=-0.5;
Defaults.WeightForAIFMeanVesses=0.3;
Options=Defaults;
MskCTC=MskC;
% [CVI, BinCVI, Bin2CVI DataToFitA]=ChooseRepVoxelsForAIFFind(Msk3D,MskE,CTC2D(:,1:NumVols),BolusStart,Options.Rep_MaxAroundBolus,Options.Rep_RatioToEnd,Options.Rep_nPerSet);
[CVI, BinCVI, Bin2CVI DataToFitA]=ChooseRepVoxelsForAIFFind(MskF,MskE,DataOut,BolusStart,Options.Rep_MaxAroundBolus,Options.Rep_RatioToEnd,Options.Rep_nPerSet);
CVI=InIdxs(CVI);
% Options=struct('SubSecRes',{[]},'MaxTDif',{[]});
%Get the relative noise (percentage out of max signal value ) of each representing voxel( rmadCTC2D is the noise of each one)
% DataNoise=rmadCTC2D(CVI);
DataNoise=EstimateNoise(DataToFitA);
% The data ( c(t) ) of the representing voxels
DataToFit=CTC2D(CVI,:);
DataToFit=DataToFitA;
% Plot2DDataOnSubfigures(11,DataToFit);

DataToFit2=DataToFit;
DataNoise=EstimateNoise(DataToFit2);
CVI2=CVI;
BinCVI2=BinCVI;
Bin2CVI2=Bin2CVI;
MskCTC2=MskCTC;
CTC2D2=CTC2D;
disp('Found representative voxels');
% Plot2DDataOnSubfigures(12,DataToFit2)
%%
ShowData=false;
if(ShowData)
    gfig(15);
    subplot(1,2,1);
    plot(DataToFit2');
    subplot(1,2,2);
    plot(NormalizeByRows(DataToFit2)');
    hold on;
    plot(mean(NormalizeByRows(DataToFit2))','k','LineWidth',2)
    plot(median(NormalizeByRows(DataToFit2))','g','LineWidth',2)
end
%%
Tmp=max(FBrainMask,[],3);
F=find(max(Tmp,[],2));
GoodRows=F(1):F(end);
F=find(max(Tmp,[],1));
GoodCols=F(1):F(end);
PrepareFN=[WorkingP 'AfterPrepare4D.mat'];
load(PrepareFN,'BadSlicesF2');
GoodSlices=setdiff(1:size(FBrainMask,3),BadSlicesF2);

Sts=0:ReqDt:(ReqDt*(NumVols-1));
GoodTs=Sts;
GoodTIdxs=1:numel(Sts);

USStr='';
RepVoxFN=[WorkingP 'RepVox' USStr '.mat'];

Mx=max(CTC2D,[],2);
SMx=sort(Mx);
MaxAmp=SMx(numel(Mx)-10);
% Show the Chosen
%     Tmp=ShowImageWithPoints(1122,CT1,MskCTC,CVI,GoodRows,GoodCols,GoodSlices);
Tmp=ShowImageWithPoints(1122,FinalT1,MskCTC2,CVI2,GoodRows,GoodCols,GoodSlices);
% Raw2Nii(Tmp*1,[WorkingP 'ChosenVoxelsForAIFFinding' USStr '.nii'],'float32', MeanFN);
% RelaxFN=[WorkingP 'Relax.mat'];
close(1122);
%% Find AIF
CVI=CVI2;
DataToFit=DataToFit2;
DataToFit2=CTC2D2(CVI,:);

DataToFit=DataToFit2;
DataNoise=EstimateNoise(DataToFit2);
% Plot2DDataOnSubfigures(12,DataToFit2)
disp('Preparing for AIF extraction');
%%
Sts=0:ReqDt:(ReqDt*(NumVols-1));
GoodTs=Sts;
GoodTIdxs=1:numel(Sts);
USStr='';
AIFFinderFN=[WorkingP 'Sim_AIFFindData' regexprep(USStr,'\.','_') '.mat'];
MaxAmp=max(DataToFit(:));
SampleTs=Sts;
TimeBetweenDCEVols=TimeBetweenDCEVolsMin*60;
alpha=1.050;beta=0.1685;s=38.078;tau=0.483;
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
nSVols=NumVols;
%     [PKOut OutAIFParam]=AIFTryf9(WorkingP,DataToFit2(:,GoodTIdxs),MaxAmp,DataNoise,ParamAIFCoeff,nSVols,TimeBetweenDCEVolsFinal,Options,false,AIFFinderFN,GoodTs,1000);
[PKOut OutAIFParam HAIF]=AIFTryf9(WorkingP,DataToFit2(:,GoodTIdxs(1:NumVols)),MaxAmp,DataNoise,ParamAIFCoeff,nSVols,TimeBetweenDCEVolsFinal,Options,false,AIFFinderFN,GoodTs(1:NumVols),1000);
close(1000:1005);
disp('Extracted the AIF');
%%
OrigCTCFN=[WorkingP 'AfterCTC' USStr '.mat'];
disp('Loading..');
OrigCTCData=load(OrigCTCFN);
UnderSampling=1;
% DCET1_PKf;
%%
Hdt=0.5/60;
TimeVec=Sts;
AIF_Parker9t=@(x,t) AIF_Parkerg3( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8),x(9),max(TimeVec))*x(2);
HSampleTs=(0:Hdt:SampleTs(end));
ExtractedHAIF=AIF_Parker9t(OutAIFParam,HSampleTs);

InspectedAIFParamsFN=[BasePath 'InspectedAIFParams.mat'];
if(exist(InspectedAIFParamsFN,'file'))
    Tmp=load(InspectedAIFParamsFN);
    OrigOutAIFParam=Tmp.InspectedParams;
else
    PKMFN   = [WorkingP 'PKM' USStr '.mat'];
    Tmp=load(PKMFN,'OutAIFParam');
    OrigOutAIFParam=Tmp.OutAIFParam;
end

OrigHAIF=AIF_Parker9t(OrigOutAIFParam,HSampleTs);
%%
MeanArt=mean(NormalizeByRows(DataToFit2))';
MedArt=median(NormalizeByRows(DataToFit2))';
gfig(101012);plot(HSampleTs,NormalizeByRows(OrigHAIF),'k','LineWidth',2);hold on;
plot(HSampleTs,NormalizeByRows(ExtractedHAIF),'b','LineWidth',2);
plot(Sts,NormalizeByRows(MeanArt'),'r','LineWidth',1);
plot(Sts,NormalizeByRows(MedArt'),'m','LineWidth',1);
title(['Corr: ' num2str(getKthElement(corrcoef(OrigHAIF,ExtractedHAIF),2))]);
legend({'True','Extracted','MeanArt','MedArt'});
disp('Finished');
%%
close(101012);
save([BasePath 'SimResults.mat'],'MeanArt','MedArt','OrigHAIF','ExtractedHAIF','OrigOutAIFParam','HSampleTs','OutAIFParam','Sts');
%% debug for start
% SimsC=load([BasePath 'SimCTCClean.mat'],'Sims');SimsC=SimsC.Sims;
% q=10000;gfig;plot(CTC2D(q,:),'*');hold on;plot(SimsC(q,:)*AIFAmpCoeff,'.r');
% T12Dx=1./(Sims*AIFAmpCoeff+repmat(1./FinalT1(MskC),1,size(Sims,2)));
% % To signal Out=SPGRfM(T1s,PDs,FAs,TR)
% Sig=zeros(size(Sims));
% for i=1:size(Sims,2)
%     disp(i);
%     Sig(:,i)=SPGRfM(T12Dx(:,i),SM0,SDCEFA,SDCETR);
% end
% 
% T12D=zeros(size(Sims));
% for i=1:size(Sims,2)
%     disp([2 i]);
%     T12D(:,i)=DCEVolToT1(FA1D,M01D,Sig(:,i));
% end
% q=122000;gfig;plot(T12D(q,:),'*');hold on;plot(T12Dx(q,:),'.r');
