BaseP='C:\STRDCE\John\Database\DCEOut\';
D=dir(BaseP);
D=D(3:end);
D=D([D.isdir]);
D={D.name}';
%%
t=0;
Tbl=zeros(numel(D),2);
for CurSub=1:numel(D)
    disp(num2str(CurSub));
    % CurSub=2;
    clearvars -regexp .* -except D CurSub BaseP t Tbl
    
    T13dFN=[BaseP D{CurSub} '\Relaxometry\Set_Num_1\T13DNFA.nii'];
    if(exist(T13dFN,'file'))
        t=t+1;
    end
    if(~exist(T13dFN,'file'))
        continue;
    end
    % end
    % %%    %
    
    
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
    
    Baseline_SegP=SPM_Segment([BaseP D{CurSub} '\Baseline.nii'],false,[],false);
    CleanBaseFN=dir([BaseP D{CurSub} '\CleanedNo*.nii']);
    CleanBaseFN=CleanBaseFN.name;
    BaselineCleaned_SegP=SPM_Segment([BaseP D{CurSub} '\' CleanBaseFN],false,[],false);
    
    BaseCSeg3D(:,:,:,1)=loadniidata([BaselineCleaned_SegP 'c1ForSeg.nii'])/256;
    BaseCSeg3D(:,:,:,2)=loadniidata([BaselineCleaned_SegP 'c2ForSeg.nii'])/256;
    BaseCSeg3D(:,:,:,3)=loadniidata([BaselineCleaned_SegP 'c3ForSeg.nii'])/256;
    [Tmp, BaseCSeg3DAll]=max(BaseCSeg3D(:,:,:,1:3),[],4);
    BaseCSeg3DAll(~EBasicMask)=0;
    
    
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
    
    CSFMask=T1CSeg3DAll==3 & T1CSeg3D(:,:,:,3)>0.99 & E2BasicMask==1;
    WMMask=T1CSeg3DAll==2 & T1CSeg3D(:,:,:,2)>0.99 & E2BasicMask==1;
    CSFMask=BaseCSeg3DAll==3 & BaseCSeg3D(:,:,:,3)>0.99 & E2BasicMask==1;
    WMMask=BaseCSeg3DAll==2 & BaseCSeg3D(:,:,:,2)>0.99 & E2BasicMask==1;
    
    T1CSeg3DAllx=T1CSeg3DAll;
    T1CSeg3DAllx=BaseCSeg3D;
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
    %%
    figure(1001);clf;
    GoodSlices=MidSli;
    GoodSlices=7;
    for i=1:numel(GoodSlices)
        CurSli=GoodSlices(i);
        I=squeeze(T1CleanedMasked(:,:,CurSli));
%         I=I*830/Tbl(CurSub,2);
        MaxV=4000;
        ClrM=[1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1];
        IRGB=repmat(min(1,I/MaxV),[1 1 3]);
        for tt=1:3
%             BW2 = bwmorph(squeeze(T1CSeg3DAll(:,:,CurSli)==tt),'remove');
            BW2 = bwmorph(squeeze(BaseCSeg3DAll(:,:,CurSli)==tt),'remove');
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
    %%
    saveas(1001,['C:\T1test\' D{CurSub} '.png']);
    saveas(1001,['C:\T1test\' D{CurSub} '.fig']);
    close(1001);
end