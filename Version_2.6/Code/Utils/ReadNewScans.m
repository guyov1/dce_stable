function AddedShortInfos = ReadNewScans(BasePath,Force)

AddedUIDs = {};
if(nargin<2)
    Force=0;
end

if (filesep ~='/') % If not Unix
    if(BasePath(1)==filesep && BasePath(2) ~=filesep)
        BasePath=[filesep BasePath];
    end
end

DD=dirdfs(BasePath);
if(sum([DD.isdir])>0)
    FDD=MultiConcat([BasePath filesep],{DD([DD.isdir]).name});
else
    FDD={BasePath};
end
for i=1:length(FDD)
    FDD{i}=regexprep(FDD{i},'\\+',filesep);
    if(FDD{i}(1)=='\' && FDD{i}(2)~='\')
        FDD{i}=['\' FDD{i}];
    end
    DCMFile{i}=getDCMFile(FDD{i});
    G(i)=~isempty(DCMFile{i});
end

%CClass=Class(G);
CDCMFile=DCMFile(G);
CPaths=FDD(G);
Infos.aaaaaaaaa=1;
    
% Infos=load(getComputerParams('InfosFN'),'Infos');
% Infos=Infos.Infos;
disp('Adding new Infos');
for i=1:length(CPaths)
    disp(['Reading ' num2str(i)]);
    Tmp=dicominfo(CDCMFile{i});
    UID=UID2FieldName(Tmp.SeriesInstanceUID);
    if(isfield(Infos,UID))
        disp(['Scan ' UID ' already exist in the database']);
        if(Force)
            disp(['Forcing ' CPaths{i}]);
            Infos(1).(UID)=Tmp;
            AddedUIDs{end+1}=Tmp.SeriesInstanceUID;
        end
    else
        Infos(1).(UID)=Tmp;
        AddedUIDs{end+1}=Tmp.SeriesInstanceUID;
    end
end
if(isfield(Infos,'aaaaaaaaa'))
    Infos=rmfield(Infos,'aaaaaaaaa');
end

ShortInfos=Infos2ShortInfos(Infos);
FNs=fieldnames(ShortInfos);
for i=1:length(FNs)
    ShortInfos.(FNs{i}).Remark=[];
    ShortInfos.(FNs{i}).Usable=1;
    ShortInfos.(FNs{i}).Path=fileparts(ShortInfos.(FNs{i}).Filename);
    ShortInfos.(FNs{i}).Class=ClassifySeq(ShortInfos.(FNs{i}).SeriesDescription,1);
    if(isfield(ShortInfos.(FNs{i}).PatientName,'FamilyName'))
        ShortInfos.(FNs{i}).Name=ShortInfos.(FNs{i}).PatientName.FamilyName;
    else
        ShortInfos.(FNs{i}).Name='None';
    end
    ShortInfos.(FNs{i}).Project=GetProjectName(ShortInfos.(FNs{i}).Path);
    if(isfield(ShortInfos.(FNs{i}),'SeriesTime') &&isfield(ShortInfos.(FNs{i}),'SeriesDate'))
        ShortInfos.(FNs{i}).SeriesDateTime=[ShortInfos.(FNs{i}).SeriesDate '_' ShortInfos.(FNs{i}).SeriesTime];
    else
        ShortInfos.(FNs{i}).SeriesDateTime='None';
    end
end

InfosC=struct2cell(ShortInfos);
AddedShortInfos=[InfosC{:}];
