% function TrgFN=CoregWrite(FN,MatFN,Force)
function TrgFN=CoregWrite(FN,MatFN,Force,TrgP,DispW,DefiningSpaceFN)
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


if(iscell(FN))
    [A S]=fileparts(FN{1});
    if(nargin<4)
        TrgP=A;
    end
    nVols=numel(FN);
    TrgFN=cell(1,nVols);
    TMat=load(MatFN,'TMat');
    TMat=TMat.TMat;
    for i=1:nVols
        if(~exist('DefiningSpaceFN','var'))
            TrgFN{i}=CoregWrite(FN{i},TMat{i},Force,TrgP,DispW);
        else
            TrgFN{i}=CoregWrite(FN{i},TMat{i},Force,TrgP,DispW,DefiningSpaceFN{i});
        end
    end
    return;
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
TrgFN=[TrgP filesep 'Coreged_' S '.nii'];
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




if(ischar(MatFN))
    TMat=load(MatFN,'TMat');
    TMat=TMat.TMat;
else
    TMat=MatFN;
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

insertTMatToNii([TempPath  '/ToReslice.nii'],TMat);
cd(TempPath);
% a=load([GetComputerParams('BaseDir')  '/Code/SPM_precofigures/Coreg_Reslice.mat']);
% a=load([getComputerParams('spmpreconfpath') 'Coreg_Reslice.mat']);
a=load(['Coreg_Reslice.mat']);
if(NoDefSpace)
    a.matlabbatch{1}.spm.spatial.coreg.write.source{1}=DefiningSpaceFN; % Target
    a.matlabbatch{1}.spm.spatial.coreg.write.ref{1}=[TempPath  '/ToResliceTrg.nii']; % Source
else
    a.matlabbatch{1}.spm.spatial.coreg.write.ref{1}=DefiningSpaceFN; % Target
    a.matlabbatch{1}.spm.spatial.coreg.write.source{1}=[TempPath  '/ToReslice.nii']; % Source
end
spm_jobman('run',a.matlabbatch);
TmpName= license('inuse');
h=findobj('Name',['SPM8 (' TmpName(1).user '): Graphics']);
close(h);
movefile([TempPath  '/rToReslice.nii'],TrgFN);
% delete([TempPath  '/Temp/Coreg/*.*']);
cd(CurP);
 [SUCCESS,MESSAGE,MESSAGEID] =  rmdir(TempPath,'s');
 if (~SUCCESS)
     display(MESSAGE)
 end
 
 
 