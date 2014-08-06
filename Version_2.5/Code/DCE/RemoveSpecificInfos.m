BPath='C:\DATA\DCE6min_control\DCEnadav\\';
APath=strrep(BPath,filesep,[filesep filesep]);
TmpInfos=getSpecificInfos('xPath',APath);
TmpFs={'shortinfosfn','infosfn'};
for jj=1:numel(TmpFs)
    a=load(getComputerParams(TmpFs{jj}));
    load(getComputerParams(TmpFs{jj}));
    FldNm=fieldnames(a);
    if(numel(FldNm)>1)
        disp('123');
    end
    for ii=1:numel(TmpInfos)
        eval([FldNm{1} '=rmfield(' FldNm{1} ',UID2FieldName(TmpInfos(ii).SeriesInstanceUID));']);
    end
    save(getComputerParams(TmpFs{jj}),FldNm{1});
end