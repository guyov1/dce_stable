BaseP='C:\DATA\AlonFriedman_Philips\AllGood\';
%%
D=dir([BaseP '*.']);
D=D(~[D.isdir]);
for i=1:numel(D)
    disp(i);
    movefile([BaseP D(i).name],[BaseP D(i).name '.dcm']);
end
%%
D=dir([BaseP '*.dcm']);
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