BaseBaseP='D:\DCE\DataBase2\';
D=dir(BaseBaseP);
D=D([D.isdir]);
D=D(3:end);
StudyNames={D.name}';
StudyNames=StudyNames(cellNumel(StudyNames)>5);
%%
for i=1:numel(StudyNames)
    disp(['Working on ' num2str([i numel(SbvtudyNames)]) ' ' StudyNames{i}]);
    clearvars -except StudyNames i BaseBaseP D
    try
        StudyName=StudyNames{i};
        WorkOnSimData;
    catch
    end
end
%%
Passed=0;
for i=1:numel(StudyNames)
    BasePath=[BaseBaseP StudyNames{i} filesep];
    Passed=Passed+(exist([BasePath 'SimResults.mat'],'file')>0);
end
%%
Cur=0;
gfig(101012);
for i=1:numel(StudyNames)numel(StudyNames)
    BasePath=[BaseBaseP StudyNames{i} filesep];
    if(~exist([BasePath 'SimResults.mat'],'file'))
        continue;
    end
    Cur=Cur+1;
    gsubplot(Passed,Cur);
    Tmp=load([BasePath 'SimResults.mat']);
    plot(Tmp.HSampleTs,NormalizeByRows(Tmp.OrigHAIF),'k','LineWidth',2);hold on;
    plot(Tmp.HSampleTs,NormalizeByRows(Tmp.ExtractedHAIF),'b','LineWidth',2);
    % plot(Sts,NormalizeByRows(MeanArt'),'r','LineWidth',1);
    % plot(Sts,NormalizeByRows(MedArt'),'m','LineWidth',1);
    title([StudyNames{i} ' Corr: ' num2str(getKthElement(corrcoef(Tmp.OrigHAIF,Tmp.ExtractedHAIF),2))]);
end
Show also MeanArt and MedArd
%% Find arteries
for CurStudyToWorkOn=1:numel(StudyNames)
    disp(['Working on ' num2str([CurStudyToWorkOn numel(StudyNames)]) ' ' StudyNames{CurStudyToWorkOn}]);
    clearvars -except StudyNames CurStudyToWorkOn BaseBaseP D
    try
        StudyName=StudyNames{CurStudyToWorkOn};
        D2=dir([BaseBaseP StudyName filesep 'DCE' filesep]);
        D2=D2([D2.isdir]);
        ShortStudyName=D2(cellNumel({D2.name})==13).name;
        BasePath=[BaseBaseP StudyName filesep 'DCE' filesep ShortStudyName filesep];
        WorkingP=BasePath;
        UnderSampling=1;

        FindArteriesFromPath;
    catch
    end
end
%%
for CurStudyToWorkOn=12:numel(StudyNames)
    disp(['Working on ' num2str([CurStudyToWorkOn numel(StudyNames)]) ' ' StudyNames{CurStudyToWorkOn}]);
    clearvars -except StudyNames CurStudyToWorkOn BaseBaseP D Datas
    StudyName=StudyNames{CurStudyToWorkOn};
    D2=dir([BaseBaseP StudyName filesep 'DCE' filesep]);
    D2=D2([D2.isdir]);
    ShortStudyName=D2(cellNumel({D2.name})==13).name;
    BasePath=[BaseBaseP StudyName filesep 'DCE' filesep ShortStudyName filesep];
    WorkingP=BasePath;
    WorkingPx=WorkingP;
    UnderSampling=1;
    USStr='';
    DebugSuffix='_ForDB';
    
    RepVoxFN=[WorkingP 'RepVox' USStr '.mat'];
    if(~exist(RepVoxFN,'file'))
        continue;
    end
    
    CTCFN=[WorkingP 'AfterCTC' USStr '.mat'];
    disp('Loading..');
    load(CTCFN,'CT1');
    load(CTCFN,'GoodRows','GoodCols','GoodSlices');
    
    WorkingP=WorkingPx;
    PrepareFN=[WorkingP 'AfterPrepare4D.mat'];
    MeanFN=[WorkingP 'DCEMean.nii'];
    
    
    a=load(RepVoxFN,'DataToFit2','CVI2','MskCTC2','CTC2D2');
    RepVoxFN=[WorkingP 'RepVox' USStr DebugSuffix '.mat'];
    b=load(RepVoxFN,'DataToFit2','CVI2','MskCTC2','CTC2D2');
    
    A1=intersect(a.CVI2,b.CVI2);
    A2=setdiff(a.CVI2,b.CVI2);
    A3=setdiff(b.CVI2,a.CVI2);
    Datas{CurStudyToWorkOn}={CT1,a.MskCTC2,GoodRows,GoodCols,GoodSlices,A1,A2,A3};
%     IAll{CurStudyToWorkOn}=ShowImageWithPointsClr(CT1,a.MskCTC2,a.CVI2,GoodRows,GoodCols,GoodSlices,TmpClr);
end
DatasM=cat(1,Datas{:});
%%
figure(827412);clf;
for CurStudyToWorkOn=1:numel(StudyNames)
    TmpClr=ones(1,32);
    IRGB3=ShowImageWithPointsClr(DatasM{CurStudyToWorkOn,1},DatasM{CurStudyToWorkOn,2},DatasM{CurStudyToWorkOn,6},DatasM{CurStudyToWorkOn,3},DatasM{CurStudyToWorkOn,4},DatasM{CurStudyToWorkOn,5},TmpClr);
    subplot(4,4,CurStudyToWorkOn);
    montage(IRGB3);
end
%%
For each of the smaller sets, find the AIF
Display: Inspected AIF, AIF from inspected RepVox, AIF from auto RepVox