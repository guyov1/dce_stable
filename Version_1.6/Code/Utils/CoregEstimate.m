% function [TrgMatFN ]=CoregEstimate(SFN,TFN,Force)
function [TrgMatFN ]=CoregEstimate(SFN,TFN,Force)
if(~exist('Force','var'))
    Force=false;
end
TmpName= license('inuse');
TmpName=TmpName(1).user;
TempPath=[getComputerParams('temppath') filesep 'Temp_CoregEst_' TmpName '_' regexprep(datestr(now),'[-: ]','_') '_' num2str(floor(rand*1000))];
mkdir(TempPath);

% delete([GetComputerParams('BaseDir') '/Temp/Coreg/*.*']);

% In unix, run the system cp command with no "-p" because it gives an
% error when the destination is in another computer so source and dest
% files have different owner
if (filesep == '/') % Unix
    system(['cp -f ' SFN ' ' [TempPath filesep 'Src.nii']]);
    system(['cp -f ' TFN ' ' [TempPath filesep 'Trg.nii']]);
    
else  % Windows
    
    copyfile(SFN,[TempPath filesep 'Src.nii'],'f');
    copyfile(TFN,[TempPath  filesep 'Trg.nii'],'f');
end


[A, S]=fileparts(SFN);
[QQQ, T]=fileparts(TFN);
CurP=pwd;
TrgMatFN=[A filesep 'CoregEstimate_' S '_to_' T '.mat'];
% TrgPSFN=[A filesep 'CoregEstimate_' S '_to_' T '.ps'];
if(exist(TrgMatFN,'file') && ~Force)
    disp(['Already exist ' TrgMatFN]);
    return;
end
cd(TempPath);
% a=load([GetComputerParams('BaseDir') '/Code/SPM_precofigures/Coreg_Estimate.mat']);
% a=load([getComputerParams('spmpreconfpath') 'Coreg_Estimate.mat']);
a=load(['Coreg_Estimate.mat']);
a.matlabbatch{1}.spm.spatial.coreg.estimate.ref{1}=[TempPath filesep 'Trg.nii']; % Target
a.matlabbatch{1}.spm.spatial.coreg.estimate.source{1}=[TempPath filesep 'Src.nii']; % Source
spm_jobman('run',a.matlabbatch);
TmpName= license('inuse');
h=findobj('Name',['SPM8 (' TmpName(1).user '): Graphics']);
close(h);
% D=dir([TempPath '/*.ps']);
% PSFN=['/local_data/Temp/Coreg/' D(1).name];
TMat=getTMatFromNii([TempPath filesep 'Src.nii']);
% delete([GetComputerParams('BaseDir') '/Temp/Coreg/*.*']);
% copyfile(PSFN,TrgPSFN,'f');
save(TrgMatFN,'TMat');
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
