%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    Set paths according to work/home computer  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get PC name
PC_name = getenv('ComputerName');

% Get Matlab user name
TmpName= license('inuse');
TmpName=TmpName(1).user;

% Remove old CompParams.mat
if (exist('CompParams.mat', 'file'))
    %path_2_remove = which 'CompParams.mat';
    delete( which('CompParams.mat') );
end


%% Set MRICron Path
% Home - Windows
if ( strcmpi(PC_name,'GUYOV-PC') )

    BaseBaseP='C:\users\';
    %mricronpath='C:\Users\Guy\Google Drive\Thesis\Matlab\Tools\mricron\';
    mricronpath='H:\Guy\Dropbox\University\Msc\Thesis\Matlab_New\Tools\mricron\';

    % Lab - Windows
% elseif ( strcmpi(PC_name,'FMRI-TOKYO') )
% 
%     BaseBaseP='C:\users\';
%     %mricronpath='C:\Users\guyn\Google Drive\Thesis\Matlab\Tools\mricron\';
%     mricronpath='C:\Users\guyn\Dropbox\Thesis\Matlab_New\Tools\mricron\';

    % Gilad's Default
    
    % Lab - Windows 2
elseif ( strcmpi(PC_name,'FMRI-HERZL') )

    BaseBaseP='D:\guyn\';
    %mricronpath='C:\Users\guyn\Google Drive\Thesis\Matlab\Tools\mricron\';
    mricronpath='D:\guyn\Dropbox\Thesis\Matlab_New\Tools\mricron\';

   % Lab - T9
elseif ( strcmpi(PC_name,'FMRI-T9') )

   %BaseBaseP='\\fmri-guy2\Dropbox\University\Msc\Thesis\Matlab_New\';
   BaseBaseP= '\\FMRI-GUY2\Dropbox\University\Msc\Thesis\SourceForge\Development\';

    %mricronpath='C:\Users\guyn\Google Drive\Thesis\Matlab_New\Tools\mricron\';
    %mricronpath = '\\fmri-guy2\Dropbox\University\Msc\Thesis\Matlab_New\Tools\mricron\';    
    mricronpath = '\\FMRI-GUY2\Dropbox\University\Msc\Thesis\SourceForge\Development\';
    
    % Gilad's Default
    
%     % Lab - T3 (Unix)
% elseif (hostname == 'fmri-t3')
% %elseif ( strcmpi(PC_name,'FMRI-T3') )
% 
%     %BaseBaseP='\\fmri-guy2\D\Dropbox\University\Msc\Thesis\';
%     %mricronpath='C:\Users\guyn\Google Drive\Thesis\Matlab_New\Tools\mricron\';
%     %mricronpath='\\fmri-herzl\D\guyn\Dropbox\University\Msc\Thesis\Matlab_New\Tools\mricron\';    
%     
%     BaseBaseP   = '/mnt/users/Moran/dcescripts/Stable_Versions/code/';
%     mricronpath = '/mnt/users/Moran/dcescripts/Stable_Versions/code/Tools/mricron/';
%     
    % Gilad's Default
    
else

    %BaseBaseP='D:\users\';
    %mricronpath='D:\users\gilad\Tools\mricron\';
    BaseBaseP=[pwd filesep];
 
end

% Add user's directory
if ( ~exist([BaseBaseP TmpName filesep 'Database'], 'file') )
    mkdir([BaseBaseP TmpName filesep 'Database']);
end

addpath(genpath([BaseBaseP 'Code'])); % this adds the DCE code 

% If using stable versions, Tools is one directory up
if ( strfind(pwd,'Version_')) 
    addpath(genpath([BaseBaseP '../' 'Tools'])); % this adds code for some utilities
else
    addpath(genpath([BaseBaseP 'Tools'])); % this adds code for some utilities
end



addpath(genpath([BaseBaseP TmpName]));

savepath([BaseBaseP TmpName filesep 'PathDefForDCE.m']);
mricronpath=[BaseBaseP 'Tools' filesep 'mricron' filesep];

%% Set Computer Parameters Struct

% Create struct for computer parameters
CompParams=struct();

% Home - Windows
if ( strcmpi(PC_name,'GUYOV-PC') )

    %save('C:\Users\Guy\Google Drive\Thesis\Matlab_New\Code\Utils\CompParams.mat','CompParams');
    save('H:\Guy\Dropbox\University\Msc\Thesis\Matlab_New\Code\Utils\CompParams.mat','CompParams');

    % Lab - Windows
% elseif ( strcmpi(PC_name,'FMRI-TOKYO') )
% 
%     %save('C:\Users\guyn\Google Drive\Thesis\Matlab_New\Code\Utils\CompParams.mat','CompParams');
%     save('C:\Users\guyn\Dropbox\Thesis\Matlab_New\Code\Utils\CompParams.mat','CompParams');

    % Lab - Windows 2
elseif ( strcmpi(PC_name,'FMRI-HERZL') )

    %save('C:\Users\guyn\Google Drive\Thesis\Matlab\Code\Utils\CompParams.mat','CompParams');
    save('D:\guyn\Dropbox\Thesis\Matlab_New\Code\Utils\CompParams.mat','CompParams');

    % Lab - T9
elseif ( strcmpi(PC_name,'FMRI-T9') )

    %save('\\fmri-guy2\Dropbox\University\Msc\Thesis\Matlab_New\Code\Utils\CompParams.mat','CompParams');
    save( '\\FMRI-GUY2\Dropbox\University\Msc\Thesis\SourceForge\Development\Code\Utils\CompParams.mat','CompParams');
    % Gilad's Default
else


    CompParams_File_Path =[BaseBaseP TmpName filesep 'Database' filesep 'CompParams.mat'];

    
    save(CompParams_File_Path,'CompParams');

end


%% Set Data's Path

% Home - Windows
if ( strcmpi(PC_name,'GUYOV-PC') )
    %BaseP= 'C:\Users\Guy\Google Drive\Thesis\Matlab\Data';
    BaseP= 'H:\Guy\Dropbox\University\Msc\Thesis\Matlab_New\Data\';

    % Lab - Windows
% elseif ( strcmpi(PC_name,'FMRI-TOKYO') )
% 
%     %BaseP= 'C:\Users\guyn\Google Drive\Thesis\Matlab\Data';
%     BaseP= 'C:\Users\guyn\Dropbox\Thesis\Matlab_New\Data\';
%     
    % Lab - Windows
elseif ( strcmpi(PC_name,'FMRI-HERZL') )

    %BaseP= 'C:\Users\guyn\Google Drive\Thesis\Matlab\Data';
    BaseP= 'D:\guyn\Dropbox\Thesis\Matlab_New\Data\';
   
    % Lab - T9
elseif ( strcmpi(PC_name,'FMRI-T9') )
 
    %BaseP= '\\fmri-guy2\Dropbox\University\Msc\Thesis\Matlab_New\Data\';
    BaseP= '\\FMRI-GUY2\Dropbox\University\Msc\Thesis\SourceForge\Development\Data\';
    
    % Gilad's Default
else
    BaseP=[BaseBaseP TmpName filesep 'Database' filesep];

end

%% Activate "setComputerParamM" function on the paths we defined

setComputerParamM('projectnamepathdepth',4);
setComputerParamM('basepath',BaseP);
setComputerParamM('temppath',[BaseP 'Temp']);
setComputerParamM('infosfn',[BaseP 'Infos' filesep 'Infos.mat']);
setComputerParamM('remarkinfosfn',[BaseP 'Infos' filesep 'RemarkInfos.mat']);
setComputerParamM('shortinfosfn',[BaseP 'Infos' filesep 'ShortInfos.mat']);
setComputerParamM('veryshortinfosfn',[BaseP 'Infos' filesep 'T9_VeryShortInfos.mat']);
[idum,hostname]= system('hostname');
F=find(hostname=='.');
if(~isempty(F))
    hostname=hostname(1:F(1)-1);
end
setComputerParamM('name',hostname);
SPMP=[fileparts(which('spm')) filesep];
TPMP=[SPMP 'tpm' filesep];
setComputerParamM('tpm',{[TPMP 'grey.nii'];[TPMP 'white.nii'];[TPMP 'csf.nii']});
setComputerParamM('nniftispaths',1);
setComputerParamM('username',TmpName);

if (~exist(getComputerParams('basepath'), 'file'))
    mkdir(getComputerParams('basepath'));
end
if (~exist([getComputerParams('basepath') 'Infos'], 'file'))
    mkdir([getComputerParams('basepath') 'Infos']);
end
if (~exist(getComputerParams('temppath'), 'file'))
    mkdir(getComputerParams('temppath'));
end

setComputerParamM('mricronpath',mricronpath);
% setComputerParamM('spmpreconfpath','D:\users\cgilad\Code\Utils\SPM_precofigures\');

display('-I- DCEInit finished successfully!');