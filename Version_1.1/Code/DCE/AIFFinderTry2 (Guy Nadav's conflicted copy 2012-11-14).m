if(Options.MakeNoBATManualArtAnalysis)
    FindAIFByManualArtNoBAT(WorkingP,false);
end

if(Options.MakeBATManualArtAnalysis)
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

    DataNoise=rmadCTC2D(CVI);
    DataToFit=ArtCTCs;
    PKMFNman=[WorkingP 'ManualArtBAT_AIF2.mat'];
    MaxAmp=max(DataToFit(:));
%     Mx=max(CTC2D,[],2);
%     SMx=sort(Mx);
%     MaxAmp=SMx(numel(Mx)-10);
%
    AIFFinderFNManual=[WorkingP 'AIFFindDataManualArt2.mat'];
    Options.EM_Num_Of_Iterations=0;
    [PKOutManual OutAIFParamManual]=AIFTryf(DataToFit,MaxAmp,DataNoise,ParamAIFCoeff,nSVols,TimeBetweenDCEVols,Options,true,AIFFinderFNManual);
    save(PKMFNman,'OutAIFParamManual','DataToFit');
end

if(~Options.MakeBATAutoArtAnalysis)
    return;
end
[CVI BinCVI Bin2CVI]=ChooseRepVoxelsForAIFFind(Idx3D>0,CTC2D,BolusStart,Options.Rep_MaxAroundBolus,Options.Rep_RatioToEnd,Options.Rep_nPerSet);

% Options=struct('SubSecRes',{[]},'MaxTDif',{[]});

DataNoise=rmadCTC2D(CVI);
DataToFit=CTC2D(CVI,:);
PKMFN=[WorkingP 'PKM' USStr '.mat'];

Mx=max(CTC2D,[],2);
SMx=sort(Mx);
MaxAmp=SMx(numel(Mx)-10);
%%
AIFFinderFN=[WorkingP 'AIFFindData' USStr '.mat'];
[PKOut OutAIFParam]=AIFTryf(DataToFit,MaxAmp,DataNoise,ParamAIFCoeff,nSVols,TimeBetweenDCEVols,Options,false,AIFFinderFN);
save(PKMFN,'OutAIFParam','DataToFit');