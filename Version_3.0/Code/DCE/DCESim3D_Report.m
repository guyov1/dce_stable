BaseP='\\fmri-t9\users\Moran\DCE\ForSim\';
D=dir(BaseP);
D=D([D.isdir]);
D=D(3:end);
D=D(~strhas({D.name},'SmVl'));
D=D(~strhas({D.name},'LiHa'));
%%
for p=1:numel(D)
    try
        BATStr='WithBAT';
        NoBATStr='WithoutBAT';
        disp('a');
        WorkingP=[BaseP D(p).name filesep];
        USStr='';
        SimFN=[WorkingP 'Sim.mat'];
        SimPKM3DFNBAT=[WorkingP 'SimPKM3D' USStr '_' BATStr '.mat'];
        SimPKM3DFNNoBAT=[WorkingP 'SimPKM3D' USStr '_' NoBATStr '.mat'];
        a=load(SimFN,'SBATs','SPKs');
        b=load(SimPKM3DFNBAT,'PKs');
        c=load(SimPKM3DFNNoBAT,'PKs');
        
        CTCFN=[WorkingP 'AfterCTC' '.mat'];
        d=load(CTCFN);
        Msk=d.DBrainMask;
        Img=d.CT1;
        
        M=[a.SPKs(:,1:4) b.PKs(:,[1 4 2 3]) c.PKs(:,[1 4 2 3])];
        I=getKrandomSamples(size(M,1),1000);
        I=1:size(M,1);
        
        D1=M(I,5:8)-M(I,1:4);
        D2=M(I,9:12)-M(I,1:4);
        DD=abs(D1)-abs(D2);
        R1=D1./M(I,1:4);
        R2=D2./M(I,1:4);
        %
        for MIdx=2:4
            Mxs=[0 1.3 0.01 0.02];
            DMxs=[0 0.1 0.01 0.0031];
            K0=Msk*0;
            K0(Msk)=M(:,MIdx);
            K1=Msk*0;
            K1(Msk)=M(:,MIdx+4);
            K2=Msk*0;
            K2(Msk)=M(:,MIdx+8);
            DK1=K1-K0;
            DK2=K2-K0;
            
            RK1=abs(DK1)./K0;
            RK2=abs(DK2)./K0;
            B=isfinite(RK1) & isfinite(RK2) & K0>0.005;
            H1=histc(RK1(B),0:0.02:2);
            H2=histc(RK2(B),0:0.02:2);
            Res(p,MIdx,:)=[mean(RK1(B)) mean(RK2(B))];
            Hists{p,MIdx}=[H1 H2];
%             Raw2Nii(K0,[WorkingP 'KOrig_' num2str(MIdx) '.nii'],'float32', [WorkingP 'DCEMean.nii']);
%             Raw2Nii(K1,[WorkingP 'KBATCorrect_' num2str(MIdx) '.nii'],'float32', [WorkingP 'DCEMean.nii']);
%             Raw2Nii(K2,[WorkingP 'KBATNoCorrect_' num2str(MIdx) '.nii'],'float32', [WorkingP 'DCEMean.nii']);
%             Raw2Nii(DK1,[WorkingP 'KDiffToBAT_' num2str(MIdx) '.nii'],'float32', [WorkingP 'DCEMean.nii']);
%             Raw2Nii(DK2,[WorkingP 'KDiffToNoBAT_' num2str(MIdx) '.nii'],'float32', [WorkingP 'DCEMean.nii']);

%                     SliI=8;
%                     figure;
%                     subplot(2,3,1);imagesc(mritransform(K0(d.GoodRows,d.GoodCols,SliI)),[0 Mxs(MIdx)]);
%                     set(gca,'XTick',[]);set(gca,'YTick',[]);
%                     subplot(2,3,2);imagesc(mritransform(K1(d.GoodRows,d.GoodCols,SliI)),[0 Mxs(MIdx)]);
%                     set(gca,'XTick',[]);set(gca,'YTick',[]);
%                     subplot(2,3,3);imagesc(mritransform(K2(d.GoodRows,d.GoodCols,SliI)),[0 Mxs(MIdx)]);
%                     set(gca,'XTick',[]);set(gca,'YTick',[]);
%                     subplot(2,3,5);imagesc(mritransform(DK1(d.GoodRows,d.GoodCols,SliI)),[-DMxs(MIdx) DMxs(MIdx)]);
%                     set(gca,'XTick',[]);set(gca,'YTick',[]);
%                     subplot(2,3,6);imagesc(mritransform(DK2(d.GoodRows,d.GoodCols,SliI)),[-DMxs(MIdx) DMxs(MIdx)]);
%                     set(gca,'XTick',[]);set(gca,'YTick',[]);
        end
        %
        
        
%         B=isfinite(M(I,2)) & M(I,2)>0.01 & isfinite(R1(:,2));
%         B2=isfinite(M(I,2)) & M(I,2)>0.01 & isfinite(R2(:,2));
%         [h,p]=ttest2(abs(R1(B,2)),abs(R2(B2,2)),'Vartype','unequal','Tail','right');
%         
%         figure;
%         nH=100;
%         subplot(1,3,1);[H X]=hist(DD(:,2),nH);bar(X,H);ax=axis;axis([ax(1:3) getKthElement(sort(H,'descend'),3)]);
%         subplot(1,3,2);[H X]=hist(DD(:,3),nH);bar(X,H);ax=axis;axis([ax(1:3) getKthElement(sort(H,'descend'),3)]);
%         subplot(1,3,3);[H X]=hist(DD(:,4),nH);bar(X,H);ax=axis;axis([ax(1:3) getKthElement(sort(H,'descend'),3)]);
        disp(['Finished writing files for ' num2str(p) ' ' BATStr ' ' D(p).name]);
    catch
        disp(['Error report in ' num2str(p) ' ' BATStr ' ' D(p).name]);
    end
end
save('Sim3DReport_MADandHists.mat','Res','Hists');
%%
for p=1:numel(D)
    for MIdx=2:4
        tmp=MIdx*2-3;
        MAD(p,tmp:tmp+1)=Res(p,MIdx,:);
    end
end