AddToLog(WorkingP,'d_00','\\subsection*{Finding the AIF}');

%% If chosen by the user, find manually the AIF
if(Options.MakeNoBATManualArtAnalysis)
    FindAIFByManualArtNoBAT(WorkingP,false);
end

if(Options.MakeBATManualArtAnalysis)
    AddToLog(WorkingP,'d_am','MakeBATManualArtAnalysis');
    ManualArt3D=loadniidata([WorkingP 'manualArt.nii']);
    F1=find(ManualArt3D);
    F2=find(Msk2);
    CVI=find(ismember(F2,F1));
    ArtCTCs=CTC2D(CVI,:);
    MeanArtCTC=mean(ArtCTCs,1);
    MeanArtCTC=MeanArtCTC*max(ArtCTCs(:))./max(MeanArtCTC);
    MeanArtCTC2=mean(NormalizeByRows(ArtCTCs),1);
    MeanArtCTC2=MeanArtCTC2*max(ArtCTCs(:))./max(MeanArtCTC2);

    figure(100);clf;subplot(2,2,1);plot(SampleTs,ArtCTCs');subplot(2,2,2);plot(SampleTs,NormalizeByRows(ArtCTCs)');
    subplot(2,2,3);plot(SampleTs,ArtCTCs');hold on;plot(SampleTs,MeanArtCTC,'k','LineWidth',2);plot(SampleTs,MeanArtCTC2,'m','LineWidth',2);
    subplot(2,2,4);plot(SampleTs,NormalizeByRows(ArtCTCs)'*max(ArtCTCs(:)));hold on;plot(SampleTs,MeanArtCTC,'k','LineWidth',2);plot(SampleTs,MeanArtCTC2,'m','LineWidth',2);
    gprint(100,[WorkingP 'manualArt.png']);
    close(100);
    AddToLog(WorkingP,'d_am1','manualArt','manualArt.png');

    DataNoise = rmadCTC2D(CVI);
    DataToFit = ArtCTCs;
    PKMFNman  = [WorkingP 'ManualArtBAT_AIFx.mat'];
    MaxAmp    = max(DataToFit(:));
%     Mx=max(CTC2D,[],2);
%     SMx=sort(Mx);
%     MaxAmp=SMx(numel(Mx)-10);
%
    AIFFinderFNManual                = [WorkingP 'AIFFindDataManualArtx.mat'];
    Options.EM_Num_Of_Iterations     = 0;
    [PKOutManual, OutAIFParamManual] = AIFTryf(WorkingP,DataToFit,MaxAmp,DataNoise,ParamAIFCoeff,nSVols,TimeBetweenDCEVols,Options,true,AIFFinderFNManual);
    save(PKMFNman,'OutAIFParamManual','DataToFit');
end

if(~Options.MakeBATAutoArtAnalysis)
    return;
end


%%
% load(RelaxFN,'GoodSlices','GoodRows','GoodCols');
Tmp=max(FBrainMask,[],3);
F=find(max(Tmp,[],2));
GoodRows=F(1):F(end);
F=find(max(Tmp,[],1));
GoodCols=F(1):F(end);
load(PrepareFN,'BadSlicesF2');
GoodSlices=setdiff(1:size(FBrainMask,3),BadSlicesF2);
%%  Get the representing voxels data and noise
[CVI, BinCVI, Bin2CVI]=ChooseRepVoxelsForAIFFind(Idx3D>0,CTC2D,BolusStart,Options.Rep_MaxAroundBolus,Options.Rep_RatioToEnd,Options.Rep_nPerSet);

% Options=struct('SubSecRes',{[]},'MaxTDif',{[]});
%Get the relative noise (percentage out of max signal value ) of each representing voxel( rmadCTC2D is the noise of each one)
DataNoise=rmadCTC2D(CVI);
% The data ( c(t) ) of the representing voxels
DataToFit=CTC2D(CVI,:);

GoodTs=Sts;
GoodTIdxs=1:numel(Sts);
if (Num_Of_T1_Maps>1)
    NumBefore=size(Additional_before_main,2);
    AfterIdxs=find(Additional_T1_Maps_Time_Diff_Sets>0)';
%     GoodTIdxs=[(NumBefore+1):(numel(ts)+NumBefore) numel(ts)+SingleAngleIdxs-1];
%     GoodTs=[Sts(1:numel(ts)) Additional_T1_Maps_Time_Diff_Sets(SingleAngleIdxs)'*86400/60];
    
    GoodTIdxs=[(NumBefore+1):(numel(ts)+NumBefore) numel(ts)+AfterIdxs-1];
    GoodTs=[Sts(1:numel(ts)) Additional_T1_Maps_Time_Diff_Sets(AfterIdxs)'*86400/60];
    % figure;plot(GoodTs,DataToFit(1:10:end,GoodTIdxs)','.-');
end

PKMFN=[WorkingP 'PKM' USStr '.mat'];

Mx=max(CTC2D,[],2);
SMx=sort(Mx);
MaxAmp=SMx(numel(Mx)-10);
% Show the Chosen
Tmp=Idx3D*0>1;
Tmp2=find(Idx3D>0)*0>1;
Raw2Nii(Tmp*1,[WorkingP 'ChosenVoxelsForAIFFinding.nii'],'float32', MeanFN);
Tmp2(CVI)=true;
Tmp(Idx3D>0)=Tmp2;
TmpD=imdilate(Tmp,strel('disk',4));
RelaxFN=[WorkingP 'Relax.mat'];

clear IRGB3
for i=1:numel(GoodSlices)
    I=squeeze(CT1(:,:,GoodSlices(i)));
    I=min(4000,I)/4000;
    Tmp2=squeeze(TmpD(:,:,GoodSlices(i)));
    I(Tmp2)=0;
    Tmp2=squeeze(Tmp(:,:,GoodSlices(i)));
    I(Tmp2)=1;
    IRGB(:,:,1)=I;
    I(Tmp2)=0;
    IRGB(:,:,2)=I;
    Tmp2=squeeze(TmpD(:,:,GoodSlices(i)));
    I(Tmp2)=1;
    Tmp2=squeeze(Tmp(:,:,GoodSlices(i)));
    I(Tmp2)=0;
    IRGB(:,:,3)=I;
    IRGB3(:,:,:,i)=IRGB(GoodRows,GoodCols,:);
end
figure(1122);clf;
montage(mritransform(IRGB3))
title('Chosen for AIF estimation');
saveas(1122,[WorkingP 'ChosenForAIFEstimation.png']);
saveas(1122,[WorkingP 'ChosenForAIFEstimation.fig']);
close(1122);
AddToLog(WorkingP,'yd_at1','Chosen for AIF estimation','ChosenForAIFEstimation.png');
%% 
AIFFinderFN=[WorkingP 'AIFFindData' regexprep(USStr,'\.','_') '.mat'];

if(Options.MakeBATPopArtAnalysis)
%     Parker, Geoff JM, et al. "Experimentally?derived functional form for a population?averaged high?temporal?resolution arterial input function for dynamic contrast?enhanced MRI." Magnetic Resonance in Medicine 56.5 (2006): 993-1000.
    A1=0.809;A2=0.330;T1=0.17046;T2=0.364;
    sig1=0.055;sig2=0.134;alpha=1.064;beta=0.166;
    s=37.772;tau=0.482;
    % Before there were really slighltly different numbers, from other place?
    % A~mM.min, T~min; sig~min, alpha~mM, beta,s~1/min, tau~min
    % Time stamp (in minutes) for every temporal point
    TimeBetweenDCEVolsMin=TimeBetweenDCEVols/60;
    SampleTs=((1:nSVols)-1)*TimeBetweenDCEVolsMin;
    MaxAmp=max(DataToFit(:));
    T1x=BolusStartSec/60;
    % Population average c(t) according to Parker's article.
    C=AIF_Parker(SampleTs,A1,sig1,T1x,A2,sig2,T2+T1x-T1,alpha,beta,s,tau+T1x-T1);
    C=(C./max(C))*MaxAmp;
    MeanVessel=C;
    AIF_Parker9t=@(x,t) AIF_Parkerg3( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8),x(9),max(SampleTs))*x(2);
    AIF_Parker9tx=@(x,t) AIF_Parker9t([x(1:8) x(9)/MaxAmp],t).*MaxAmp;
    
    MinFirstBolusSig=Options.MinFirstBolusStd; %2; % seconds
    LB=[0.1 0     1/60 0.1     0   0.01  0.01 0]';
    UB=[10  1.5   0.25                10      1   3    2   0.3]';
    % X0=[T1  A1  sig1                T2Delta A2  sig2 Alpha Beta]';
    X0=[1   1     0.15                1 0.3 0.3   0.4 0.2]';
    CostFuncMeanVessel=@(x) gCost(MeanVessel,AIF_Parker9tx(min(max(x,[LB; 0]'),[UB; 1]'),SampleTs),'RMS');
%     OldParams1=fminsearch(CostFuncMeanVessel,[X0' MeanVessel(end)]);
    X1=[T1x   1     0.05                T2 0.03 0.1   0.1 0.3 MeanVessel(end)];
    OldParams1=fminsearch(CostFuncMeanVessel,X1);
    OldParams1=min(max(OldParams1,[LB; 0]'),[UB; 1]');
%%
    TmpAIF=AIF_Parker9tx(OldParams1,SampleTs);
    TmpAIF1=AIF_Parker9tx(X1,SampleTs);
%     figure(101);clf;plot(SampleTs,C,'k','LineWidth',3);hold on;
%     plot(SampleTs,TmpAIF,'r*');
    OutAIFParam=OldParams1;
%%
else
    [PKOut OutAIFParam]=AIFTryf9(WorkingP,DataToFit(:,GoodTIdxs),MaxAmp,DataNoise,ParamAIFCoeff,nSVols,TimeBetweenDCEVolsFinal,Options,false,AIFFinderFN,GoodTs);
end

% Save the stage .mat file
save(PKMFN,'OutAIFParam','DataToFit');