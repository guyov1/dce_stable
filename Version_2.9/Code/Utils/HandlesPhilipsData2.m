BaseP='\\fmri-t9\users\df-1_backup\General_Group_data\DCE\Case from Verona\DI_GIOR\REQ2013\DCE_Verona_Philips\';
%%
DD=dir([BaseP '*.']);
DD=DD(3:end);
for j=1:numel(DD)
    CurP=[BaseP DD(j).name filesep]
    D=dir([CurP '*.']);
    D=D(~[D.isdir]);
    for i=1:numel(D)
        disp(i);
        movefile([CurP D(i).name],[CurP D(i).name '.dcm']);
    end
end
%%
D=dirdfs(BaseP);
DNames={D.name};
D=D(strhas(DNames,'dcm'));
for i=1:numel(D)
    disp(i);
    Infos{i}=dicominfo([BaseP D(i).name]);
end
save([BaseP 'Infos.mat'],'Infos');
% InfosS=[Infos{:}];
%%
for i=1:numel(D)
    UID=UID2FieldName(Infos{i}.SeriesInstanceUID);
    AllInfos.(UID)=Infos{i};
end
ShortInfos=Infos2ShortInfos(AllInfos);
ShortInfosA=struct2array(ShortInfos);
%%
gDicom2Nifti(BaseP,[BaseP 'Out\A.nii']);
%%
load([BaseP 'Infos.mat']);
%%
s=warning('off','MATLAB:MKDIR:DirectoryExists');
%%
for i=1:numel(Infos)
    if(~(strhas(Infos{i}.SeriesDescription,'uptake') || strhas(Infos{i}.SeriesDescription,'DESPOT')))
        continue;
    end
    disp(i);
    TrgP=[BaseP 'Out2' filesep Infos{i}.SeriesDescription filesep];
    mkdir(TrgP);
    if(Infos{i}.FileSize<40000)
        continue;
    end
    copyfile(Infos{i}.Filename,[TrgP 'MR0' num2str(Infos{i}.InstanceNumber) '.dcm']);
end