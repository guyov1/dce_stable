function UpdateShortInfos(Force)
if(~exist('Force','var'))
    Force=0;
end
tempDirName = dir(getComputerParams('ShortInfosFN'));
% <<<<<<< .mine
% <<<<<<< .mine
% if(isfield(tempDirName,'datenum'))
%     LastUpdate = tempDirName.datenum;
% else
%     LastUpdate = floor( datenum(tempDirName.date) );
% end
% CurDate = floor( datenum(date) );
% if(CurDate-LastUpdate==0)
%     disp('Already updated');
%     return;
% =======
% if(~isempty(tempDirName) && ~Force)
%     LastUpdate = floor( datenum(tempDirName.date) );
%     CurDate = floor( datenum(date) );
%     if(CurDate-LastUpdate==0)
%         disp('Already updated');
%         return;
%     end
% % >>>>>>> .r567
% end
% =======
% if(isfield(tempDirName,'datenum'))
%     LastUpdate = tempDirName.datenum;
% else
%     LastUpdate = floor( datenum(tempDirName.date) );
% end
% CurDate = floor( datenum(date) );
% if(CurDate-LastUpdate==0 && ~Force)
%     disp('Already updated');
%     return;
% end
% >>>>>>> .r617
disp('Updating short infos');
Q=load(getComputerParams('InfosFN'));
ShortInfos=Infos2ShortInfos(Q.Infos);
FNs=fieldnames(ShortInfos);
if(exist(getComputerParams('remarkinfosfn'),'file'))
    load(getComputerParams('remarkinfosfn'));
end
if(~exist('Remarks','var'))
    Usable=struct([]);
end
if(~exist('Remarks','var'))
    Remarks=struct([]);
end
for i=1:length(FNs)
    ShortInfos.(FNs{i}).Remark=[];
    ShortInfos.(FNs{i}).Usable=1;
    if(isfield(Usable,FNs{i}))
        ShortInfos.(FNs{i}).Usable=Usable.(FNs{i});
    end
    if(isfield(Remarks,FNs{i}))
        ShortInfos.(FNs{i}).Remark=Remarks.(FNs{i});
    end
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
save(getComputerParams('ShortInfosFN'),'ShortInfos');
UpdateVeryShortInfos;