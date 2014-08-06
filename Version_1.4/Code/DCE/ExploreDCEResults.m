%%
MeanVol=loadniidata(MeanFN);

load([WorkingP 'CTCMsk' USStr '.mat'],'CTC2DBigGood','MskCTCGood');
load(CTCFN,'DBrainMask');

MskCTCGood3D=DBrainMask;
MskCTCGood3D(MskCTCGood3D)=MskCTCGood;

Idx3DBig=NaN(size(Msk));
Idx3DBig(MskCTCGood3D)=1:sumn(MskCTCGood3D);

TTls={'BATfinal' 'VpFinal' 'KtransFinal' 'Kepfinal' 'VeFinal' 'RSSFinal' 'RSS0' 'RSS1' 'RSS2' 'RSS3' 'F1v0' 'F2v1' 'F3v2' 'BAT1' 'Vp1' 'BAT2' 'Vp2' 'Ktrans2' 'BAT3' 'Vp3' 'Ktrans3' 'Kep3' 'Ve3' 'WhichModel' 'TProb1' 'TProb2' 'TProb3' 'TProbFinal'};

TimeBetweenDCEVolsMin=TimeBetweenDCEVols/60;
InterpolationFactor=ceil(TimeBetweenDCEVols);
SampleTs=((1:nSVols)-1)*TimeBetweenDCEVolsMin;

HInterpolationFactor=ceil(InterpolationFactor*2);
Hdt=TimeBetweenDCEVolsMin/HInterpolationFactor;
HSampleTs=0:Hdt:SampleTs(end);

ThreeSec=ceil(Options.MaxTDif_ForWholeVOI/(Hdt*60));
TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;
nTDif=numel(TDif);
AIF_Parker8t=@(x,t) AIF_Parkerg2( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8))*x(2);

if(Options.MakeBATManualArtAnalysis)
    PKM3DFNman=[WorkingP 'PKM3DManualx.mat'];
    load(PKM3DFNman);
    HInterpolationFactor=ceil(InterpolationFactor*Options.SubSecondResolution);
    Hdt=TimeBetweenDCEVolsMin/HInterpolationFactor;
    HSampleTs=0:Hdt:SampleTs(end);

    ThreeSec=ceil(Options.MaxTDif_ForWholeVOI/(Hdt*60));
    TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;
    nTDif=numel(TDif);

    HAIF=AIF_Parker8t(OutAIFParam,HSampleTs);

    AllVols=cell(1,29);
    for i=1:numel(TTls)
        AllVols{i+1}=MskCTCGood3D*0;
        AllVols{i+1}(MskCTCGood3D)=PKs(:,i);
    end
    AllVols{1}=MeanVol;
    AllVols{2}(AllVols{2}==0)=NaN;
    AllVols{2}(~isnan(AllVols{2}))=-60*TDif(AllVols{2}(~isnan(AllVols{2})));
    AllVols{6}(isnan(AllVols{6}) & AllVols{25}==3)=Inf;
%     AllVols{5}(isnan(AllVols{6}) & AllVols{25}==1)=0;
%     AllVols{25}=AllVols{25}+1;
    Titles=['MeanVol' TTls];

    DCEResultsGUI(AllVols,Titles,MskCTCGood3D,CTC2DBigGood,Idx3DBig,1,PKs,1,1,1,SampleTs,HSampleTs,HAIF,PKM3DFNman,TDif);

%     DCEResultsGUI({MeanVol Kep3D BAT3D Vp3D Ktrans3D Ve3D -log10(F3D) -log10(F03D) TP3D MskCTCGood3D MaxConcentration3D Noise3D (MaxRatioToBaseline-1)*100 RMS3D rRMS3D RMStoNoise3D CT1 B1},Titles,MskCTCGood3D,CTC2DBigGood,Idx3DBig,1,PKs,1,1,1,SampleTs,HSampleTs,HAIF,PKM3DFNman);

    return;
end

%%
PKM3DFN=[WorkingP 'PKM3D' USStr '.mat'];
load(PKM3DFN,'OutAIFParam','PKs');
disp(PKM3DFN);
disp(OutAIFParam);

HAIF=AIF_Parker8t(OutAIFParam,HSampleTs);

%%
Kep3D=MskCTCGood3D*0;
Kep3D(MskCTCGood3D)=PKs(:,2);
BAT3D=MskCTCGood3D*0;

BATVals=PKs(:,1);
BATVals(BATVals==0 | isnan(BATVals))=nTDif+1;
TDifx=[TDif NaN];
BAT3D(MskCTCGood3D)=-60*TDifx(BATVals);

Vp3D=MskCTCGood3D*0;
Vp3D(MskCTCGood3D)=PKs(:,3);
Ktrans3D=MskCTCGood3D*0;
Ktrans3D(MskCTCGood3D)=PKs(:,4);
Ve3D=Ktrans3D./Kep3D;
Ve3D(~MskCTCGood3D)=0;
RMS3D=MskCTCGood3D*0;
RMS3D(MskCTCGood3D)=PKs(:,6)*1000;
F3D=MskCTCGood3D*0;
F3D(MskCTCGood3D)=max(0,PKs(:,7));
F03D=MskCTCGood3D*0;
F03D(MskCTCGood3D)=max(0,PKs(:,8));
RMSs=PKs(:,6)';
TP3D=MskCTCGood3D*0;
TP3D(MskCTCGood3D)=max(0,PKs(:,9));
TP3D(TP3D==2)=-0.5;
% 
load(CTCFN,'CT1','B1');

PrepareFN=[WorkingP 'AfterPrepare4D.mat'];
load(PrepareFN,'MaxRatioToBaseline');

MaxConcentration3D=NaN(size(Msk));
MaxConcentration3D(MskCTCGood3D)=max(CTC2DBigGood,[],2)*1000;

rRMS3D=RMS3D*NaN;
rRMS3D(MskCTCGood3D)=RMSs./max(CTC2DBigGood,[],2)';

Noise2D=EstimateNoise(CTC2DBigGood);
Noise3D=RMS3D*NaN;
Noise3D(MskCTCGood3D)=Noise2D;

RMStoNoise3D=RMS3D./(MaxConcentration3D.*max(Noise3D,1e-1));


%PKs(:,1)=-60*TDif(PKs(:,1));
PKs(:,1)=-60*TDifx(BATVals);
Titles={'MeanVol' 'Kep','BAT','Vp','Ktrans','Ve' 'F' 'F0' 'TimeP' 'Msk' 'maxC' 'Noise' 'Enhancement%' 'RMS' 'rRMS3D' 'RMStoNoise' 'T1','B1'};
TTL=PKM3DFN;
DCEResultsGUI({MeanVol Kep3D BAT3D Vp3D Ktrans3D Ve3D -log10(F3D) -log10(F03D) TP3D MskCTCGood3D MaxConcentration3D Noise3D (MaxRatioToBaseline-1)*100 RMS3D rRMS3D RMStoNoise3D CT1 B1},Titles,MskCTCGood3D,CTC2DBigGood,Idx3DBig,1,PKs,1,1,1,SampleTs,HSampleTs,HAIF,TTL,TDif);
% DCEResultsGUI({MeanVol Kep3D BAT3D Vp3D Ktrans3D Ve3D F3D MskCTCGood3D
% MaxConcentration3D Noise3D (MaxRatioToBaseline-1)*100 RMS3D rRMS3D RMStoNoise3D CT1 B1},Titles,MskCTCGood3D,CTC2DBigGood,Idx3DBig,1,PKs,HSAIF,HSHConvd,Keps1I,SampleTs,HSampleTs,HAIF);