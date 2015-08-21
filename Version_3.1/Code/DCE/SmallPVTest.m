BaseP='C:\STRDCE\John\Database\DCEOut\';
D=dir(BaseP);
D=D(3:end);
D=D([D.isdir]);
D={D.name}';
%%
TODO:
V Automatic recognition of when Baseline already has CA!!
V Draw a figure of the different regions
V Extract the main ventricle median T1 and normalize it so that generally WM is around 830        
V Extract a map of WM too (for AIF normalization) and save the figure
V Extract CTC using both the new FA and B1
V Save a map of the representatives for AIF optimization
PDF report of all the process
%%
for CurSub=1:numel(D)
    a=load([BaseP D{CurSub} '\Relaxometry\Set_Num_1\NFARes.mat']);
    ONFAsC{CurSub}=a.ONFAs;
    LL(CurSub)=ONFAsC{CurSub}(end);
    if(any(a.ONFAs==0.5) || sum((a.ONFAs-floor(a.ONFAs))==0)>1)
        disp(['Problem with ' D{CurSub}]);
        Tmp=loadniidata([BaseP D{CurSub} '\DCEMainCoreged\Coreged_vol_0001.nii']);
        figure(CurSub);
        imagesc(mritransform(Tmp(:,:,floor(size(Tmp,3)/2))));
        colormap gray;
        title(D{CurSub});
    end
end
%%
% for CurSub=1:numel(D)
    CurSub=3;
    disp(num2str(CurSub));
    
    T13dFN=[BaseP D{CurSub} '\Relaxometry\Set_Num_1\T13DNFA.nii'];
    if(~exist(T13dFN,'file'))
        continue;
    end
    
    a=load([BaseP D{CurSub} '\Relaxometry\Set_Num_1\NFARes.mat']);
    ONFAs=a.ONFAs;
    
    PrepareFN=[BaseP D{CurSub} '\AfterPrepare4D.mat'];
    T1NFA_SegP=SPM_Segment(T13dFN,false,[],false);
    
    T1Seg3D(:,:,:,1)=loadniidata([T1NFA_SegP 'c1ForSeg.nii'])/256;
    T1Seg3D(:,:,:,2)=loadniidata([T1NFA_SegP 'c2ForSeg.nii'])/256;
    T1Seg3D(:,:,:,3)=loadniidata([T1NFA_SegP 'c3ForSeg.nii'])/256;
    % T1Seg3D(:,:,:,4)=loadniidata(FMaskFN);
    [Tmp, T1Seg3DAll]=max(T1Seg3D(:,:,:,1:3),[],4);
    % T1Seg3DAll(~T1Seg3D(:,:,:,4))=0;
    
    T1Cleaned=loadniidata([T1NFA_SegP 'mForSeg.nii']);
    load(PrepareFN,'BadSlicesF2');
    BasicMask=sum(T1Seg3D,4)>0.99;
    BasicMask(:,:,BadSlicesF2)=false;
    BasicMask=bwfillHoles3Dby2D(BasicMask);
    se=strel('disk',4,8);
    EBasicMask=imerode(BasicMask,se);
    EBasicMaskFN=[BaseP D{CurSub} '\Relaxometry\Set_Num_1\EBasicMask.nii'];
    Raw2Nii(EBasicMask,EBasicMaskFN,'float32', T13dFN);
    T1Seg3DAll(~EBasicMask)=0;
    
    T1NFAC_SegP=SPM_Segment([T1NFA_SegP 'mForSeg.nii'],false,[],false,EBasicMaskFN);
    
    T1CSeg3D(:,:,:,1)=loadniidata([T1NFAC_SegP 'c1ForSeg.nii'])/256;
    T1CSeg3D(:,:,:,2)=loadniidata([T1NFAC_SegP 'c2ForSeg.nii'])/256;
    T1CSeg3D(:,:,:,3)=loadniidata([T1NFAC_SegP 'c3ForSeg.nii'])/256;
    [Tmp, T1CSeg3DAll]=max(T1CSeg3D(:,:,:,1:3),[],4);
    T1CSeg3DAll(~EBasicMask)=0;
    
    MidSli=floor(size(BasicMask,3)/2);
    GoodSlices=setdiff(1:size(BasicMask,3),BadSlicesF2);
    
    T1CleanedMasked=T1Cleaned;
    T1CleanedMasked(~EBasicMask)=0;
    % end
    % %%
    E2BasicMask=imerode(BasicMask,strel('disk',16,8));
    
    RMS3DFN=[BaseP D{CurSub} '\Relaxometry\Set_Num_1\RMS3DNFA.nii'];
    RMS3D=loadniidata(RMS3DFN);
    
%     I=squeeze(T1Cleaned(:,:,MidSli));
    I=squeeze(RMS3D(:,:,MidSli));
    P=[121 154]; % I(P(1),257-P(2))
    DD=dir([BaseP D{CurSub} '\Relaxometry\Set_Num_1\Coreged_*']);
    bla=@(x) x{1};
    Tmp=regexp({DD.name}','.*_(.*)\.nii','tokens');
    SerNum=str2double(cellfun(bla,Tmp));
    [X Y]=sort(SerNum);
    for i=1:numel(Y)
        FA4DF(:,:,:,Y(i))=loadniidata([BaseP D{CurSub} '\Relaxometry\Set_Num_1\' DD(i).name]);
    end
    FA4DF(:,:,:,end+1)=loadniidata([BaseP D{CurSub} '\Baseline.nii']);
%     load(PrepareFN,'CurMainDCEInfo'); % CurMainDCEInfo.RepetitionTime
    TRs=[a.TRsF ];
    %%
    CurSli=7;
    I=squeeze(RMS3D(:,:,CurSli));
    figure(121212);imagesc(mritransform(I),[0 100]);
%     figure;imagesc(mritransform(I),[0 4000]);
    P=floor(ginput(1));
    close(121212);
    
    MaskI=E2BasicMask(:,:,CurSli);
    [FI FJ]=find(MaskI);
    for ff=1:numel(FI)
        disp(ff);
        %     Data=squeeze(FA4DF(P(1),257-P(2),MidSli,:));
        Data=squeeze(FA4DF(FI(ff),FJ(ff),CurSli,:));
        RegMat=@(T11,T12) SPGRfM(min(10000,max(0,[T11 T12])),[1 1],ONFAs,TRs);
        CostByT1s=@(T1s) gCost(Data',max(Data'/(RegMat(T1s(1),T1s(2))),0)*RegMat(T1s(1),T1s(2)),'RMS');
        [BestX, BestVal]=fminsearch(CostByT1s,[3000 800]);
        BestX=min(10000,max(0,BestX));
        PDs=max(Data'/(RegMat(BestX(1),BestX(2))),0);
        T1s2D(ff,:)=BestX;
        PDs2D(ff,:)=PDs;
    end
    %%
    TmpI=MaskI*0;
    TmpI(MaskI)=min(T1s2D,[],2);
    figure;imagesc(mritransform(TmpI),[0 5000]);
    TmpI=MaskI*0;
    TmpI(MaskI)=max(T1s2D,[],2);
    figure;imagesc(mritransform(TmpI),[0 5000])
    
    TmpI=MaskI*0;
    TmpI(MaskI)=max(PDs2D,[],2)./min(PDs2D,[],2);
    figure;imagesc(mritransform(TmpI),[0 5])
    
    [MM MI]=max(PDs2D,[],2);
    TmpI=MaskI*0;
    TmpI(MaskI)=T1s2D(sub2ind(size(PDs2D),1:size(PDs2D,1),MI'));
    figure;imagesc(mritransform(TmpI),[0 5000])
    %%
    A=SPGRfM(BestX,PDs,ONFAs,TRs);
    [T1 PD]=CalcT1byFAfw2(Data,ONFAs,TRs);
    Sim1=SPGRfM(T1,PD,ONFAs,TRs);
    figure(3030);clf;
    plot(Data,'k-');hold on;
    plot(A(1,:),'r-');
    plot(A(2,:),'m-');
    plot(sum(A),'g-');
    plot(Sim1,'b-');
    title(num2str([BestX T1],'% 2.0f'));
    xlabel(num2str([PDs PD],'% 2.0f'));
    %%
    WarningStatus=warning('off','MATLAB:rankDeficientMatrix');
    FFA2DData=Reshape4d22d(FA4DF,BasicMask>-1);
    [RelaxCG{1:3}]=CalcT1byFAfw2(FFA2DData',ONFAs,TRs);
    SimF=SPGRfM(RelaxCG{1}',RelaxCG{2}',ONFAs,TRs);
    
    [RelaxCG{1:3}]=CalcT1byFAfw2(FFA2DData',a.FAsF,TRs);
    SimO=SPGRfM(RelaxCG{1}',RelaxCG{2}',a.FAsF,TRs);
    %%
    figure(93893);clf;
    CurSli=7;
    for rr=1:size(FA4DF,4)
        TmpI=squeeze(FA4DF(:,:,CurSli,rr));
        Tmp3D=squeeze(FA4DF(:,:,:,1));
        Tmp3D(:)=SimO(:,rr);
        TmpI2=squeeze(Tmp3D(:,:,CurSli));
        gsubplot(2,size(FA4DF,4),1,rr);
        imagesc(mritransform(TmpI-TmpI2),[-50 50]);
        
        Tmp3D(:)=SimF(:,rr);
        TmpI2=squeeze(Tmp3D(:,:,CurSli));
        gsubplot(2,size(FA4DF,4),2,rr);
        imagesc(mritransform(TmpI-TmpI2),[-50 50]);
    end
    %%
    
    
    CSFMask=T1CSeg3DAll==3 & T1CSeg3D(:,:,:,3)>0.99 & E2BasicMask==1;
    WMMask=T1CSeg3DAll==2 & T1CSeg3D(:,:,:,2)>0.99 & E2BasicMask==1;
    
    T1CSeg3DAllx=T1CSeg3DAll;
    T1CSeg3DAllx(CSFMask)=4;
    T1CSeg3DAllx(WMMask)=5;
    
    Tmp=T1CleanedMasked(CSFMask);
    Tmp=Tmp(isfinite(Tmp));
    Tbl(CurSub,1)=median(Tmp(Tmp>0));
    Tmp=T1CleanedMasked(WMMask);
    Tmp=Tmp(isfinite(Tmp));
    Tbl(CurSub,2)=median(Tmp(Tmp>0));
    % end
    %
    if(Tbl(CurSub,1)>10000)
        continue;
    end
    figure(1001);clf;
    GoodSlices=MidSli;
    for i=1:numel(GoodSlices)
        CurSli=GoodSlices(i);
        I=squeeze(T1CleanedMasked(:,:,CurSli));
        MaxV=4000;
        ClrM=[1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1];
        IRGB=repmat(min(1,I/MaxV),[1 1 3]);
        for tt=1:3
            BW2 = bwmorph(squeeze(T1CSeg3DAll(:,:,CurSli)==tt),'remove');
            for kk=1:3
                TmpI=IRGB(:,:,kk);
                TmpI(BW2)=ClrM(tt,kk);
                IRGB(:,:,kk)=TmpI;
            end
        end
        for tt=4:5
            BW2 = bwmorph(squeeze(T1CSeg3DAllx(:,:,CurSli)==tt),'remove');
            for kk=1:3
                TmpI=IRGB(:,:,kk);
                TmpI(BW2)=ClrM(tt,kk);
                IRGB(:,:,kk)=TmpI;
            end
        end
        gsubplot(numel(GoodSlices),i);
        imagesc(mritransform(IRGB))
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        title(CurSli);
        xlabel(CurSub);
        ylabel(Tbl(CurSub,:));
    end
    saveas(1001,['C:\T1test\' D{CurSub} '.png']);
    saveas(1001,['C:\T1test\' D{CurSub} '.fig']);
    close(1001);
end