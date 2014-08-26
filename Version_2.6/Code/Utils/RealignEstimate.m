% function [TrgMatFN]=RealignEstimate(SFN,Force)
function [TrgMatFN]=RealignEstimate(SFN,Force,DispW)
if(~exist('Force','var'))
    Force=false;
end
if(~exist('DispW','var'))
    DispW=true;
end
[A, S]=fileparts(SFN{1});
CurP=pwd;
TrgMatFN=[A filesep 'RealignEstimate_' S '.mat'];
% TrgPSFN=[A filesep 'CoregEstimate_' S '_to_' T '.ps'];
if(exist(TrgMatFN,'file') && ~Force)
    if(DispW)
        disp(['Already exist ' TrgMatFN]);
    end
    return;
end

TmpName= license('inuse');
TmpName=TmpName(1).user;
TempPath=[getComputerParams('temppath') filesep 'Temp_RealignEst_' TmpName '_' regexprep(datestr(now),'[-: ]','_') '_' num2str(floor(rand*1000)) filesep];
mkdir(TempPath);
% TempPath=[getComputerParams('temppath') filesep 'Coreg' filesep];
% delete([TempPath '*.*']);
nVols=numel(SFN);
AFN=cell(1,nVols);
for i=1:nVols
    AFN{i}=[TempPath 'vol_' sprintf('%02.0f',i) '.nii'];
    
    % In unix, run the system cp command with no "-p" because it gives an
    % error when the destination is in another computer so source and dest
    % files have different owner
    if (filesep == '/') % Unix
        system(['cp -f ' SFN{i} ' ' AFN{i}]);
    else  % Windows
        copyfile(SFN{i},AFN{i},'f');
    end
end
cd(TempPath);
% a=load('/data/home/gilad/Desktop/Code/SPM_precofigures/Realign_Estimate.mat');
% a=load([getComputerParams('spmpreconfpath') 'Realign_Estimate.mat']);
a=load(['Realign_Estimate.mat']);
a.matlabbatch{1}.spm.spatial.realign.estimate.data{1}=AFN; % Target
spm_jobman('run',a.matlabbatch);
TmpName= license('inuse');
h=findobj('Name',['SPM8 (' TmpName(1).user '): Graphics']);
close(h);
TxtFile=[TempPath 'rp_vol_01.txt'];
Movement=importdata(TxtFile);
TMat=cell(1,nVols);
for i=1:nVols
    TMat{i}=getTMatFromNii(AFN{i});
end
% delete([TempPath '*.*']);
save(TrgMatFN,'TMat','Movement');
cd(CurP);

pause(0.5);
% Try to remove the directory 5 times before failing
for i=1:5
    [success, message, msg_id] = rmdir(TempPath,'s');
    if success
        break;
    else 
        pause(1);
        continue;
    end
end

if (~success)
    display(['Error! Unable to remove the following dir: ' TempPath]);
    display(['Error message: ' message]);
    display(['Error message ID: ' msg_id]);
    exit;
end
