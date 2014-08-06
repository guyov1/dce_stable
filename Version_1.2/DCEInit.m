%%    Set needed paths

% Enable parallel processing
%matlabpool;

% Get PC name
PC_name = getenv('ComputerName');

% Get Matlab user name
TmpName= license('inuse');
TmpName=TmpName(1).user;

% Remove old CompParams.mat
if (exist('CompParams.mat', 'file'))
    delete( which('CompParams.mat') );
end

% Set base path
BaseBaseP=[pwd filesep];

% Add the DCE code path
addpath(genpath([BaseBaseP 'Code']));

% Add user's Database directory
if ( ~exist([BaseBaseP TmpName filesep 'Database'], 'file') )
    mkdir([BaseBaseP TmpName filesep 'Database']);
end

% Adding the used utilities path
if ( strfind(pwd,'Version_')) %If using "stable versions", Tools is one directory up
    addpath(genpath([BaseBaseP '../' 'Tools']));
    % Set MRI-Cron path
    mricronpath=[BaseBaseP '../' 'Tools' filesep 'mricron' filesep];
else
    addpath(genpath([BaseBaseP 'Tools']));
    % Set MRI-Cron path
    mricronpath=[BaseBaseP 'Tools' filesep 'mricron' filesep];
end

% Add user's directory path
addpath(genpath([BaseBaseP TmpName]));

% Save the path list in user's directory
savepath([BaseBaseP TmpName filesep 'PathDefForDCE.m']);



% Set Computer Parameters Struct
CompParams=struct();
CompParams_File_Path =[BaseBaseP TmpName filesep 'Database' filesep 'CompParams.mat'];
save(CompParams_File_Path,'CompParams');

% Set Base Path as Data's Path
BaseP=[BaseBaseP TmpName filesep 'Database' filesep];

%%  Set Computer Parameters according to paths defined above

setComputerParamM('projectnamepathdepth',4);
setComputerParamM('basepath',BaseP);
setComputerParamM('temppath',[BaseP 'Temp']);
setComputerParamM('infosfn',[BaseP 'Infos' filesep 'Infos.mat']);
setComputerParamM('remarkinfosfn',[BaseP 'Infos' filesep 'RemarkInfos.mat']);
setComputerParamM('shortinfosfn',[BaseP 'Infos' filesep 'ShortInfos.mat']);
setComputerParamM('veryshortinfosfn',[BaseP 'Infos' filesep 'T9_VeryShortInfos.mat']);
setComputerParamM('mricronpath',mricronpath);
setComputerParamM('nniftispaths',1);
setComputerParamM('username',TmpName);
% setComputerParamM('spmpreconfpath','D:\users\cgilad\Code\Utils\SPM_precofigures\');

% % Get the host name
% [idum,hostname]= system('hostname');
% F=find(hostname=='.');
% if(~isempty(F))
%     hostname=hostname(1:F(1)-1);
% end
% setComputerParamM('name',hostname);

% Set SPM's needed path
SPMP=[fileparts(which('spm')) filesep];
TPMP=[SPMP 'tpm' filesep];
setComputerParamM('tpm',{[TPMP 'grey.nii'];[TPMP 'white.nii'];[TPMP 'csf.nii']});

% Create the base,temp and info directories (if does not exist already)
if (~exist(getComputerParams('basepath'), 'file'))
    mkdir(getComputerParams('basepath'));
end
if (~exist([getComputerParams('basepath') 'Infos'], 'file'))
    mkdir([getComputerParams('basepath') 'Infos']);
end
if (~exist(getComputerParams('temppath'), 'file'))
    mkdir(getComputerParams('temppath'));
end

% Display success message if finished
display('-I- DCEInit finished successfully!');