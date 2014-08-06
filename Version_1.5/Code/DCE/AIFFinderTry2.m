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

    DataNoise=rmadCTC2D(CVI);
    DataToFit=ArtCTCs;
    PKMFNman=[WorkingP 'ManualArtBAT_AIFx.mat'];
    MaxAmp=max(DataToFit(:));
%     Mx=max(CTC2D,[],2);
%     SMx=sort(Mx);
%     MaxAmp=SMx(numel(Mx)-10);
%
    AIFFinderFNManual=[WorkingP 'AIFFindDataManualArtx.mat'];
    Options.EM_Num_Of_Iterations=0;
    [PKOutManual OutAIFParamManual]=AIFTryf(WorkingP,DataToFit,MaxAmp,DataNoise,ParamAIFCoeff,nSVols,TimeBetweenDCEVols,Options,true,AIFFinderFNManual);
    save(PKMFNman,'OutAIFParamManual','DataToFit');
end

if(~Options.MakeBATAutoArtAnalysis)
    return;
end


%%  Get the representing voxels data and noise
[CVI BinCVI Bin2CVI]=ChooseRepVoxelsForAIFFind(Idx3D>0,CTC2D,BolusStart,Options.Rep_MaxAroundBolus,Options.Rep_RatioToEnd,Options.Rep_nPerSet);

% Options=struct('SubSecRes',{[]},'MaxTDif',{[]});
%Get the relative noise (percentage out of max signal value ) of each representing voxel( rmadCTC2D is the noise of each one)
DataNoise=rmadCTC2D(CVI);
% The data ( c(t) ) of the representing voxels
DataToFit=CTC2D(CVI,:);

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
load(RelaxFN,'GoodSlices','GoodRows','GoodCols');
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

AIFFinderFN=[WorkingP 'AIFFindData' USStr '.mat'];

% Get the optimal paramters (Vp, Ktrans, Kep, delta T)
[PKOut OutAIFParam]=AIFTryf(WorkingP,DataToFit,MaxAmp,DataNoise,ParamAIFCoeff,nSVols,TimeBetweenDCEVolsFinal,Options,false,AIFFinderFN);

% Save the stage .mat file
save(PKMFN,'OutAIFParam','DataToFit');