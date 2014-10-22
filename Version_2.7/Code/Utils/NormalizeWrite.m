% function TrgFN=NormalizeWrite(FN,MatFN,Force,TrgP,DispW,DefiningSpaceFN)
function TrgFN=NormalizeWrite(FN,MatFN,Force,TrgP,DispW,DefiningSpaceFN)
if(~exist('Force','var'))
    Force=[];
end
if(isempty(Force))
    Force=false;
end
if(~exist('DispW','var'))
    DispW=[];
end
if(isempty(Force))
    DispW=true;
end


TmpName= license('inuse');
TmpName=TmpName(1).user;
TempPath=[getComputerParams('temppath') filesep 'Temp_CoregWrite_' TmpName '_' regexprep(datestr(now),'[-: ]','_') '_' num2str(floor(rand*1000))];

if(~exist('DefiningSpaceFN','var'))
    DefiningSpaceFN=[TempPath  '/ToReslice.nii'];
end
NoDefSpace=strcmp([TempPath  '/ToReslice.nii'],DefiningSpaceFN);
[A S]=fileparts(FN);
if(nargin<4)
    TrgP=A;
end
if(isempty(TrgP))
    TrgP=A;
end
CurP=pwd;
FileSepStr=filesep;
if(TrgP(end)==filesep)
    FileSepStr='';
end
TrgFN=[TrgP FileSepStr 'Conormalized_' S '.nii'];
if(exist(TrgFN,'file') && ~Force)
    if(DispW)
        disp(['Already exist ' TrgFN]);
    end
    return;
end

mkdir(TempPath);

% delete([GetComputerParams('BaseDir')  '/Temp/Coreg/*.*']);

% In unix, run the system cp command with no "-p" because it gives an
% error when the destination is in another computer so source and dest
% files have different owner
if (filesep == '/') % Unix
    system(['cp -f ' FN ' ' [TempPath  '/ToReslice.nii']]);
else  % Windows
    copyfile(FN,[TempPath  '/ToReslice.nii'],'f');
end

if(NoDefSpace)
    % In unix, run the system cp command with no "-p" because it gives an
    % error when the destination is in another computer so source and dest
    % files have different owner
    if (filesep == '/') % Unix
        system(['cp -f ' FN ' ' [TempPath  '/ToResliceTrg.nii']]);
    else  % Windows
        copyfile(FN,[TempPath  '/ToResliceTrg.nii'],'f');
    end
end

cd(TempPath);
% a=load([GetComputerParams('BaseDir')  '/Code/SPM_precofigures/Coreg_Reslice.mat']);
% a=load([getComputerParams('spmpreconfpath') 'Coreg_Reslice.mat']);
a=load(['Normalize_Write.mat']);
if(NoDefSpace)
    a.matlabbatch{1}.spm.spatial.normalise.write.subj.matname{1}=MatFN; % Target
    a.matlabbatch{1}.spm.spatial.normalise.write.subj.resample{1}=[TempPath  '/ToResliceTrg.nii']; % Source
else
%     a.matlabbatch{1}.spm.spatial.normalise.write.ref{1}=DefiningSpaceFN; % Target
%     a.matlabbatch{1}.spm.spatial.normalise.write.source{1}=[TempPath  '/ToReslice.nii']; % Source
end
spm_jobman('run',a.matlabbatch);
TmpName= license('inuse');
h=findobj('Name',['SPM8 (' TmpName(1).user '): Graphics']);
close(h);
movefile([TempPath  '/wToResliceTrg.nii'],TrgFN);
% delete([TempPath  '/Temp/Coreg/*.*']);
cd(CurP);
[SUCCESS,MESSAGE,MESSAGEID] =  rmdir(TempPath,'s');
if (~SUCCESS)
    display(MESSAGE)
end