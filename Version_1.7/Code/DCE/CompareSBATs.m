BaseP='C:\DCE\John\Database\DCEOut\';
D=dir(BaseP);
D=D([D.isdir]);
D=D(3:end);
D=D(~strhas({D.name},'SmVl'));
D=D(~strhas({D.name},'LiHa'));
figure(3000);clf;
MaxP=numel(D);
for p=1:MaxP %numel(D)
    try
        WorkingP=[BaseP D(p).name filesep];
        disp(WorkingP);
        %
        USStr='';
        CTCFN=[WorkingP 'AfterCTC' '.mat'];
        PKMFN=[WorkingP 'PKM' USStr '.mat'];
        PKM3DFN=[WorkingP 'PKM3D' USStr '.mat'];
        
        SimFN=[WorkingP 'Sim.mat'];
        SimPKMFN=[WorkingP 'SimPKM' USStr '.mat'];
        SimPKM3DFN=[WorkingP 'SimPKM3D' USStr '.mat'];
        
        a=load(CTCFN);
        b=load(PKMFN,'OutAIFParam');
        c=load(a.PrepareFN,'TimeBetweenDCEVols','BolusStart');
        d=load(PKM3DFN,'PKs');
        
        load(SimFN);
        load(SimPKMFN);
        load(SimPKM3DFN);
        %
        nSVols=size(a.CTC2D,2);
        TimeBetweenDCEVolsMin=c.TimeBetweenDCEVols/60;
        InterpolationFactor=ceil(c.TimeBetweenDCEVols);
        SampleTs=((1:nSVols)-1)*TimeBetweenDCEVolsMin;
        SampleTsC{p}=SampleTs;
        
        HInterpolationFactor=ceil(InterpolationFactor*2);
        Hdt=TimeBetweenDCEVolsMin/HInterpolationFactor;
        HSampleTs=0:Hdt:SampleTs(end);
        HSampleTsC{p}=HSampleTs;
        
        AIF_Parker8t=@(x,t) AIF_Parkerg2( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8))*x(2);
        EAIF=AIF_Parker8t(OutAIFParam,HSampleTs);
        EAIFC{p}=EAIF;
        HAIFC{p}=HAIF;
        OutAIFParamC{p}=OutAIFParam;
        NewParamsC{p}=NewParams;
        PKsC{p}=PKs;
        SPKsC{p}=SPKs;
        
        % [b.OutAIFParam; NewParams; OutAIFParam]
        gsubplot(MaxP*2,p*2-1);
        plot(HSampleTs,HAIF,'k',HSampleTs,EAIF,'g',HSampleTs,EAIF*max(HAIF)/max(EAIF),'m');
        
        ThreeSec=ceil(3/(Hdt*60))*2;
        TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;
        nTDif=numel(TDif);
        
        TBATs=TDif(PKs(:,1));
        %
        
        TMsk=TBATs>min(TBATs) & TBATs<max(TBATs);
        TMsk=TMsk & PKs(:,9)'<0.95;
        Edg=unique(TBATs(TMsk));
        ToCompare=[SBATs(TMsk) TBATs(TMsk)'];
        ToCompareC{p}={SBATs TBATs PKs SPKs};
        HH=hist3(ToCompare,{Edg,Edg});
        gsubplot(MaxP*2,p*2);
        imagesc(HH);
        title(getKthElement(corrcoef(ToCompare),2));
    catch
        disp('Problem');
    end
end
disp('--a--');
%
Ps=0.5:0.01:0.99;
clear CCP MM
for p=1:MaxP-1
    try
        TBATs=ToCompareC{p}{2};
        SBATs=ToCompareC{p}{1};
        TMska=TBATs>min(TBATs) & TBATs<max(TBATs) & SBATs'>min(SBATs') & SBATs'<max(SBATs');
        for i=1:numel(Ps)
            TMsk=TMska & ToCompareC{p}{3}(:,9)'<Ps(i);
            ToCompare=[SBATs(TMsk) TBATs(TMsk)'];
            CCP(p,i)=getKthElement(corrcoef(ToCompare),2);
        end
        TMsk=TMska & ToCompareC{p}{3}(:,9)'<0.95;
        MM(p,:)=[std(SBATs(TMsk)) std(TBATs(TMsk)) CCP(p,end-4) sum(TMsk)/100000];
    catch
    end
end
GoodS=MM(:,1)>2/60;
FGoodS=find(GoodS);
nGood=numel(FGoodS);
disp('--b--');
figure(5000);clf;
plot(Ps,CCP(GoodS,:)');
%% AIFs
clear AIFCorrs
figure(1);clf;hold on;
set(gca,'FontSize',20);
for i=1:nGood
    CurP=FGoodS(i);
%     gsubplot(nGood,i);
    GraphDelay=(i-1)*2;
    plot(HSampleTsC{CurP}+GraphDelay,HAIFC{CurP}./max(HAIFC{CurP}),'k-','LineWidth',2);
    plot(HSampleTsC{CurP}+GraphDelay,EAIFC{CurP}/max(EAIFC{CurP}),'-','Color',ones(1,3)*0.5,'LineWidth',2);
    AIFCorrs(i)=gcorr(HAIFC{CurP},EAIFC{CurP});
%     title(AIFCorrs(i));
end
mAIFCorrs=mean(AIFCorrs);
title(['Mean AIF Corr: ' num2str(mAIFCorrs)]);
set(gca,'YTick',[]);
xlabel('Time (min)');
axis([0 18 0 1]);
FigFN='SimAIFs';
saveas(1,[FigFN '.fig']);
saveas(1,[FigFN '.tif']);
%% BAT
clear BATCorrs
AllSBATs=[];
AllEBATs=[];
figure(20);clf;subplot(1,2,1);
BATPThresh=0.95;
for i=1:nGood
    CurP=FGoodS(i);
    SBATs=ToCompareC{CurP}{1};
    TBATs=ToCompareC{CurP}{2}';
%     TMsk=TBATs>min(TBATs) & TBATs<max(TBATs);
%         TMsk=TMsk & PKs(:,9)'<0.95;
%     TMska=TBATs>min(TBATs) & TBATs<max(TBATs) & PKsC{CurP}(:,9)<0.95; % & SBATs>min(SBATs) & SBATs<max(SBATs)
    TMska=PKsC{CurP}(:,9)>0 & PKsC{CurP}(:,9)<BATPThresh;
    DataS=SBATs(TMska)-median(SBATs(TMska));
    DataE=TBATs(TMska)-median(TBATs(TMska));
    AllSBATs=[AllSBATs; DataS];
    AllEBATs=[AllEBATs; DataE];
    [C P]=corrcoef(DataS,DataE);
    BATCorrs(i,:)=[C(2) P(2)];
end
mBATCorrs=mean(BATCorrs(:,1));
Edg=unique(AllTBATs);
ToCompare=[AllSBATs AllEBATs];
HH=hist3(ToCompare,{Edg,Edg});
HHA=HH*0;
for i=1:nGood
    CurP=FGoodS(i);
    SBATs=ToCompareC{CurP}{1};
    TBATs=ToCompareC{CurP}{2}';
%     TMska=TBATs>min(TBATs) & TBATs<max(TBATs) & SBATs>min(SBATs) & SBATs<max(SBATs) & PKsC{CurP}(:,9)<0.95;
    TMska=TBATs>min(TBATs) & PKsC{CurP}(:,9)>0 & PKsC{CurP}(:,9)<BATPThresh; % & SBATs>min(SBATs) & SBATs<max(SBATs)
    TMskb=PKsC{CurP}(:,9)>0;
    PInMsk(i)=sum(TMska)./sum(TMskb);
%     TMska=PKsC{CurP}(:,9)<0.95;
    DataS=SBATs(TMska)-median(SBATs(TMska));
    DataE=TBATs(TMska)-median(TBATs(TMska));
    ToCompare=[DataS DataE];
    HHC{i}=hist3(ToCompare,{Edg,Edg});
    HHA=HHA+HHC{i}./sumn(HHC{i});
end
imagesc(HHA);
title(['Mean BAT Corr: ' num2str(mBATCorrs)]);
% title([mBATCorrs]);
% BAT Ps
Ps=0.9:0.01:1;
subplot(1,2,2);cla;hold on;
clear CCPX PPX
for j=1:nGood
    CurP=FGoodS(j);
    SBATs=ToCompareC{CurP}{1};
    TBATs=ToCompareC{CurP}{2}';
    TMska=TBATs>min(TBATs) & TBATs<max(TBATs) & SBATs>min(SBATs) & SBATs<max(SBATs) & ToCompareC{CurP}{3}(:,9)>0;
    for i=1:numel(Ps)
        TMsk=TMska & ToCompareC{CurP}{3}(:,9)<Ps(i);
        ToCompare=[SBATs(TMsk) TBATs(TMsk)];
        CCPX(j,i)=getKthElement(corrcoef(ToCompare),2);
        PPX(j,i)=sum(TMsk)./sum(TMska);
    end
    plot(PPX(j,:),CCPX(j,:));
end
title('BAT correlation as function of the percent of voxels included');
% plot(Ps,CCPX(GoodS,:)');

%% Vp
RandSSize=10000;
% BAT Kep Vp Ktrans Ve RMS fVal f0Val
figure(22);clf; subplot(1,4,1);hold on;
CurPK=3;
CLRsM=jet(nGood);
for i=1:nGood
    CurP=FGoodS(i);
%     gsubplot(nGood,i);
    GoodV=PKsC{CurP}(:,8)<0.1;
    nGoodV=numel(GoodV);
    DataS=SPKsC{CurP}(GoodV,CurPK);
    DataE=PKsC{CurP}(GoodV,CurPK);
    RI=getKrandomSamples(numel(DataS),min(numel(DataS),RandSSize));
%     plot(DataS(RI),DataE(RI),'.','Color',CLRsM(i,:));
    scatter(DataS(RI),DataE(RI),1,CLRsM(i,:),'filled')
    VPCorrs(i)=gCost(DataS,DataE,'Corr');
    VpTTLs{i}=(['Study# ' num2str(i) ', n=' num2str(nGoodV) ', r=' num2str(VPCorrs(i))]);
end
% legend(VpTTLs);
axis equal
MaxVp=1;
axis([0 MaxVp 0 MaxVp]);
mVPCorrs=mean(VPCorrs);
title(['Mean Vp Corr: ' num2str(mVPCorrs)]);
% Ktrans
% BAT Kep Vp Ktrans Ve RMS fVal f0Val
subplot(1,4,2);cla; hold on;
CurPK=4;
for i=1:nGood
    CurP=FGoodS(i);
%     gsubplot(nGood,i);
    MaxKtrans=2.2;
    GoodV=PKsC{CurP}(:,7)<0.1 & SPKsC{CurP}(:,CurPK)<MaxKtrans;
    nGoodV=numel(GoodV);
    DataS=SPKsC{CurP}(GoodV,CurPK);
    DataE=PKsC{CurP}(GoodV,CurPK);
    RI=getKrandomSamples(numel(DataS),min(numel(DataS),RandSSize));
%     plot(DataS(RI),DataE(RI),'.','Color',CLRsM(i,:));
    scatter(DataS(RI),DataE(RI),1,CLRsM(i,:),'filled')
    KtransCorrs(i)=gCost(DataS,DataE,'Corr');
%     title([nGoodV KtransCorrs(i)]);
    KtransTTLs{i}=(['Study# ' num2str(i) ', n=' num2str(nGoodV) ', r=' num2str(KtransCorrs(i))]);
end
% legend(KtransTTLs);
axis equal
axis([0 MaxKtrans 0 MaxKtrans]);
mKtransCorrs=mean(KtransCorrs);
title(['Mean Ktrans Corr: ' num2str(mKtransCorrs)]);
% Kep
% BAT Kep Vp Ktrans Ve RMS fVal f0Val
subplot(1,4,3);cla;hold on;
CurPK=2;
for i=1:nGood
    CurP=FGoodS(i);
%     gsubplot(nGood,i);
    GoodV=PKsC{CurP}(:,7)<0.1;
    nGoodV=numel(GoodV);
    DataS=SPKsC{CurP}(GoodV,CurPK);
    DataE=PKsC{CurP}(GoodV,CurPK);
    RI=getKrandomSamples(numel(DataS),min(numel(DataS),RandSSize));
%     plot(DataS(RI),DataE(RI),'.','Color',CLRsM(i,:));
    scatter(DataS(RI),DataE(RI),1,CLRsM(i,:),'filled')
    KepCorrs(i)=gCost(DataS,DataE,'Corr');
%     title([nGoodV KepCorrs(i)]);
    KepTTLs{i}=(['Study# ' num2str(i) ', n=' num2str(nGoodV) ', r=' num2str(KepCorrs(i))]);
end
% legend(KepTTLs);
axis equal
axis([0 5 0 5]);
mKepCorrs=mean(KepCorrs);
title(['Mean Kep Corr: ' num2str(mKepCorrs)]);
% Ve
% BAT Kep Vp Ktrans Ve RMS fVal f0Val
subplot(1,4,4);cla;hold on;
CurPK=5;
for i=1:nGood
    CurP=FGoodS(i);
%     gsubplot(nGood,i);
    GoodV=PKsC{CurP}(:,7)<0.1 & SPKsC{CurP}(:,2)>0 & SPKsC{CurP}(:,4)<MaxKtrans & SPKsC{CurP}(:,4)./SPKsC{CurP}(:,2)<10; % & SPKsC{CurP}(:,5)<8 & PKsC{CurP}(:,5)<9;
    nGoodV=numel(GoodV);
%     DataS=SPKsC{CurP}(GoodV,CurPK);
    DataS=SPKsC{CurP}(GoodV,4)./SPKsC{CurP}(GoodV,2);
    DataE=PKsC{CurP}(GoodV,CurPK);
    RI=getKrandomSamples(numel(DataS),min(numel(DataS),RandSSize));
%     plot(DataS(RI),DataE(RI),'.','Color',CLRsM(i,:));
    scatter(DataS(RI),DataE(RI),1,CLRsM(i,:),'filled')
    VeCorrs(i)=gCost(DataS,DataE,'Corr');
%     title([nGoodV VeCorrs(i)]);
    VeTTLs{i}=(['Study# ' num2str(i) ', n=' num2str(nGoodV) ', r=' num2str(VeCorrs(i))]);
end
% legend(VeTTLs);
axis equal
MaxVe=10;
axis([0 10 0 10]);
mVeCorrs=mean(VeCorrs);
title(['Mean Ve Corr: ' num2str(mVeCorrs)]);
%%





















HH=hist3([BA(FMsk) BB(FMsk)],{Edg,Edg});
figure;imagesc(HH);

FA=loadniidata('C:\DCE\John\Database\DCEOut\SmVl_20120930\F.nii');
FB=loadniidata('C:\DCE\John\Database\DCEOut\SmVl_20120930\F_3.nii');
RMS=loadniidata('C:\DCE\John\Database\DCEOut\SmVl_20120930\RMS.nii');
FThresh=10^-3;
FMsk=FA<FThresh & FB<FThresh & RMS>0;
figure;imagesc(mritransform(FMsk(:,:,2)));
%%
BA=loadniidata('C:\DCE\John\Database\DCEOut\SmVl_20120930\BAT.nii');
BB=loadniidata('C:\DCE\John\Database\DCEOut\SmVl_20120930\BAT_3.nii');
%%
FMsk=PKs(:,7)<10^-13;
FMsk=FMsk & PKs(:,1)>min(PKs(:,1)) & PKs(:,1)<max(PKs(:,1)) & EBATs'~=0;
Edg=unique(EBATs(FMsk));
ToCompare=[SBATs(FMsk) EBATs(FMsk)'];
HH=hist3(ToCompare,{Edg,Edg});
figure;imagesc(HH);
title(getKthElement(corrcoef(ToCompare),2));
%%
DB=BA-BB-median(BA(FMsk))+median(BB(FMsk));
DB(~FMsk)=NaN;
Vals=abs(DB(FMsk));
[Ns Bins]=hist(Vals,1000);
figure;plot(Bins,cumsum(Ns)./sum(Ns));
title(median(Vals));
%%
Idxs=getKrandomSamples(1:N,500);
clear TSims HTSims;
HConvd2=DCECostFuncgrT1ForConv(HAIF',SPKs(Idxs,2),HSampleTs,HConvIdxMTriB,HTriB);
for j=1:numel(Idxs)
    SAIF=interp1(HSampleTs,HAIF,SampleTs+SBATs(Idxs(j)),[],'extrap');
    HSAIF=interp1(HSampleTs,HAIF,HSampleTs+SBATs(Idxs(j)),[],'extrap');
    SHConvd2=interp1(HSampleTs,HConvd2(j,:)',SampleTs+SBATs(Idxs(j)),[],'extrap')';
    HSHConvd2=interp1(HSampleTs,HConvd2(j,:)',HSampleTs+SBATs(Idxs(j)),[],'extrap')';
    Regressors=[SAIF; SHConvd2'];
    TSims(j,:)=((Regressors')*(SPKs(Idxs(j),[3 4])'));
    Regressors=[HSAIF; HSHConvd2'];
    HTSims(j,:)=((Regressors')*(SPKs(Idxs(j),[3 4])'));
end

CHAIF=cumtrapz(HSampleTs,EAIF);
MSAIF=zeros([nTDif numel(SampleTs)]);
MCSAIF=zeros([nTDif numel(SampleTs)]);
for i=1:nTDif
    MSAIF(i,:)=interp1(HSampleTs,EAIF,SampleTs+TDif(i),[],'extrap');
    MCSAIF(i,:)=interp1(HSampleTs,CHAIF,SampleTs+TDif(i),[],'extrap');
end

TPKs = FindPKBATgAIFMuraseF_TProb(TSims,MSAIF,SampleTs,MCSAIF);
TBATs=TDif(TPKs(:,1));

clear ETSims;
for j=1:numel(Idxs)
    ESAIF=interp1(HSampleTs,HAIF,HSampleTs+TBATs(j),[],'extrap');
    ESHConvd2=interp1(HSampleTs,HConvd2(j,:)',HSampleTs+TBATs(j),[],'extrap')';
    Regressors=[ESAIF; ESHConvd2'];
    ETSims(j,:)=((Regressors')*(TPKs(j,[3 4])'));
end

TMsk=TBATs>min(TBATs) & TBATs<max(TBATs);
TMsk=TMsk & TPKs(:,9)'<0.1;
Edg=unique(TBATs(TMsk));
ToCompare=[SBATs(Idxs(TMsk)) TBATs(TMsk)'];
HH=hist3(ToCompare,{Edg,Edg});
figure;imagesc(HH);
title(getKthElement(corrcoef(ToCompare),2));

F=find(sign(TBATs').*sign(SBATs(Idxs))==-1);

i=15;
figure;plot(SampleTs,TSims(F(i),:),'k*',HSampleTs,HTSims(F(i),:),'k-',HSampleTs,ETSims(F(i),:),'b-')
%%
TTPKs = FindPKBATgAIFMuraseF_TProb(TSims(F(i),:),MSAIF,SampleTs,MCSAIF);

The noise distributes like N-dimensional gaussian
given that we have perfect fit.
P(Data|ModelA)=exp(-RSS/sqr_sig)
P(Data|ModelB)=exp(-RSS/sqr_sig)
P(Data|Model)~P(Model|Data)
The likelihood for a given
P(RSS|Perfect fit)=

P(