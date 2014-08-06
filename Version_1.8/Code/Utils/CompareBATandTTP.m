BaseP='C:\STRDCE\John\Database\DCEOut2\Compare\';
PerfusionPath{1}='C:\DATA\x\Sivi_Movement\St11_Se52_DSC-Perfusion_15secdelay_72rep';
PerFile2FN='C:\STRDCE\John\Database\DCEOut2\Compare\Perf2\GoSh_DSC_4D.nii';
A=loadniidata(PerFile2FN);
PerfPath{1}=[BaseP 'Perf1\'];
PerfPath{2}=[BaseP 'Perf2\'];
mkdir(PerfPath{1});
% gDicom2Nifti(PerfusionPath{1},[PerfPath{1} 'Perf.nii']);
PerfFile{1}=[PerfPath{1} 'Perf_0002.nii'];
PerfFile{2}=[PerfPath{2} 'Perf_0002.nii'];
% Raw2Nii(squeeze(A(:,:,:,2)),PerfFile{2},'float32', PerFile2FN);
OutP{1}='C:\STRDCE\John\Database\DCEOut\SiAv_20130213\';
OutP{2}='C:\STRDCE\John\Database\DCEOut2\GoSh_20121231\';
for i=1:2
    BaselineFN{i}=[OutP{i} 'Baseline.nii'];
    CoregFN{i}=CoregEstimate(PerfFile{i},BaselineFN{i},false);
end
TTPFN{1}='C:\DATA\forBATev\SIRI_AVNER_Analysis\dsc\xdsc_ttplr.nii';
TTPFN{2}='C:\DATA\forBATev\GOLAN_SHOSHANA_Analysis\dsc\xdsc_ttplr.nii';
% TTPFN{1}='C:\DATA\forBATev\SIRI_AVNER_Analysis\dsc\xdsc_oMTTlr.nii';
cPerfP{1}=[BaseP 'cPerf1\'];
cPerfP{2}=[BaseP 'cPerf2\'];
%%
for i=1:2
    cTTPFN{i}=RemoveDoubleFilesep(CoregWrite(TTPFN{i},CoregFN{i},false,cPerfP{i},false,BaselineFN{i}));
    BAT{i}=loadniidata([OutP{i} 'AutoArtBAT\BATfinal.nii']);
    BATP{i}=loadniidata([OutP{i} 'AutoArtBAT\TProbFinal.nii']);
    TTP{i}=round(flipdim(loadniidata(cTTPFN{i}),1));
    load([OutP{i} 'AfterCTC.mat'],'GoodSlices','EBrainMask','BadSlicesF2','GoodRows','GoodCols');
    GoodSlicesC{i}=GoodSlices;
    BadSlicesF2C{i}=BadSlicesF2;
    GoodRowsC{i}=GoodRows;
    GoodColsC{i}=GoodCols;
    se=strel('disk',4,8);
    EEBrainMask=imerode(EBrainMask,se);
    EEEBrainMask=imerode(EEBrainMask,se);
    EEEEBrainMask=imerode(EEEBrainMask,se);
    EEEEEBrainMask=imerode(EEEEBrainMask,se);
    % Mask{1}=loadniidata([OutP{1} 'BrainMask.nii']);
    Mask{i}=EEEBrainMask;
    Mask{i}(:,:,BadSlicesF2C{i})=false;
    Mask2{i}=EEEEEBrainMask;
    Mask2{i}(:,:,BadSlicesF2C{i})=false;
end
%%
% h = fspecial('disk', 9);
h=fspecial('gaussian',19,14);
Filter=@(x) imfilter(x,h);
NFilter=@(x) FilterWithNaNs(x,Filter,0.5);

for i=1:2
    BATS{i}=Apply2DFuncBySlices(BAT{i},NFilter);
%     BATS{i}=BAT{i};
%     for s=1:size(TTP{i},3)
%         BATS{i}(:,:,s)=imfilter(BAT{i}(:,:,s), h);
%     end
end
%
Names={'SiAv','GoSh'};
figure(3000);clf;
for i=1:2
    MaskA{i}=Mask2{i} & isfinite(BATS{i}) & BATS{i}>0 & TTP{i}>6 & TTP{i}<25 & BATP{i}>0;%  & BAT{i}<=15;% & BAT{i}>5;
    Vals{i}=[TTP{i}(MaskA{i}) round(BATS{i}(MaskA{i}))];
    Edges{i,1}=unique(Vals{i}(:,1));
    Edges{i,2}=unique(Vals{i}(:,2));
    H{i}=hist3(Vals{i},'Edges',Edges(i,:));
    gsubplot(2,i);
    imagesc(H{i});
    [C P]=corrcoef(Vals{i});
    title([Names{i} ', Corr: ' num2str(-C(1,2)) ' num voxels:' num2str(sumn(MaskA{i}))]);
    xlabel('BAT');
    ylabel('TTP');
end
% MaximizeSaveCloseAndAddToLog(BaseP,'Corrs');
%%
for i=1:2
    figure(10+i);clf;
    M=numel(GoodSlicesC{i});
    % M=1;
    nRows=ceil(numel(GoodSlicesC{i})*2/4);
    for ii=1:M
        CurSli=GoodSlicesC{i}(ii);
        subplot(nRows,4,ii*2-1);imagesc(mritransform((-BATS{i}(GoodRowsC{i},GoodColsC{i},CurSli)+max(BATS{i}(:))).*squeeze(MaskA{i}(GoodRowsC{i},GoodColsC{i},CurSli))));
        set(gca,'XTick',[],'YTick',[]);
        if(ii<=2) title('BAT'); end
        ylabel(ii);
        subplot(nRows,4,ii*2);imagesc(mritransform(TTP{i}(GoodRowsC{i},GoodColsC{i},CurSli).*squeeze(MaskA{i}(GoodRowsC{i},GoodColsC{i},CurSli))),[0 20]);
        set(gca,'XTick',[],'YTick',[]);
        if(ii<=2) title('TTP'); end
    end
    set(gcf,'Position',figposition([0 0 100 100]));
%     MaximizeSaveCloseAndAddToLog(BaseP,Names{i});
end
%% Output for dafna
for i=1:2
    Raw2Nii(BATS{i},['C:\STRDCE\BATTTP\BATS' num2str(i) '.nii'],'float32', cTTPFN{i});
    Raw2Nii(BAT{i},['C:\STRDCE\BATTTP\BAT' num2str(i) '.nii'],'float32', cTTPFN{i});
    Raw2Nii(TTP{i},['C:\STRDCE\BATTTP\TTP' num2str(i) '.nii'],'float32', cTTPFN{i});
end