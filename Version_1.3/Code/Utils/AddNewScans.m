function AddedShortInfos = AddNewScans(BasePath,Force)

AddedUIDs = {};
if(nargin<2)
    Force=0;
end

% BasePath='C:\Users\giora\Contacts\Desktop\U\MRI\Phantom4';
DD=dirdfs(BasePath);
if(sum([DD.isdir])>0)
    FDD=MultiConcat([BasePath filesep],{DD([DD.isdir]).name});
else
    FDD={BasePath};
end
for i=1:length(FDD)
    DCMFile{i}=getDCMFile(FDD{i});
    G(i)=~isempty(DCMFile{i});
end
%CClass=Class(G);
CDCMFile=DCMFile(G);
CPaths=FDD(G);
disp('Loading Infos');
if(exist(getComputerParams('InfosFN'),'file'))
    load(getComputerParams('InfosFN'),'Infos');
else
    Infos.aaaaaaaaa=1;
end

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
save(getComputerParams('InfosFN'),'Infos');
% <<<<<<< .mine
% UpdateShortInfos;
%
% if(numel(AddedUIDs)==0)
%     disp('Nothing added');
%     return;
% end
% temp1 = MultiConcat('|',AddedUIDs);
% temp1 = strcat(temp1{:});
% temp1 = [ '(zzz' temp1 ')' ];
% =======
UpdateShortInfos(1);
if(isempty(AddedUIDs))
    temp1='zzz';
    disp('Nothing added');
    AddedShortInfos = getSpecificInfos( 'xSeriesInstanceUID', temp1 );
    return;
else
    temp1 = MultiConcat('|',AddedUIDs);
    temp1 = strcat(temp1{:});
    temp1 = [ '(zzz' temp1 ')' ];
end
% >>>>>>> .r571
AddedShortInfos = getSpecificInfos( 'xSeriesInstanceUID', temp1 );

disp('Finished');
disp('...');
disp('...');
warning( 'off', 'MATLAB:nonIntegerTruncatedInConversionToChar' );
A1 = mat2cell( num2str((1:length(AddedShortInfos))'),ones(length(AddedShortInfos),1));
if(~isempty(AddedShortInfos))
    if(any(gIsEmpty({AddedShortInfos.SeriesNumber})))
        AddedShortInfos(gIsEmpty({AddedShortInfos.SeriesNumber})).SeriesNumber=-777;
    end
end
A2 = mat2cell(num2str( [ AddedShortInfos.SeriesNumber]'),ones(length(AddedShortInfos),1));
A3 = {AddedShortInfos.SeriesDescription}';
disp( 'Info No. -- Series No. -- Series Description' );
disp( strcat( A1, ' --', A2 , ' -- ', A3 ) )