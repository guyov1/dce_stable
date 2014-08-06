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
