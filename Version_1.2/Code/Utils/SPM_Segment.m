% function [TrgP]=SPM_Segment(SFN,Force,TrgP,DispW, MaskFN)
function [TrgP]=SPM_Segment(SFN,Force,TrgP, DispW, MaskFN)
if(~exist('Force','var'))
    Force=[];
end
if(isempty(Force))
    Force=false;
end
if(~exist('DispW','var'))
    DispW=true;
end
[A S]=fileparts(SFN);
if(~exist('TrgP','var'))
    TrgP=[];
end
if(isempty(TrgP))
    TrgP=[A filesep 'Seg_' S filesep];
end
CurP=pwd;
if(exist(TrgP,'dir') && ~Force)
    if(DispW)
        disp(['Already exist ' TrgP]);
    end
    return;
end
disp(['SPM Segmentation start for ' SFN ' ' datestr(now)]);

TmpName= license('inuse');
TmpName=TmpName(1).user;
TempPath=[getComputerParams('temppath') '_SPMSegment_' TmpName '_' regexprep(datestr(now),'[-: ]','_') '_' num2str(floor(rand*1000))];
mkdir(TempPath);

try
    % delete([GetComputerParams('BaseDir')   '/Temp/Coreg/*.*']);
    AFN=[TempPath   filesep 'ForSeg.nii'];
    copyfile(SFN,AFN,'f');
    cd(TempPath);
%     cd([GetComputerParams('BaseDir')   '/Temp/Coreg']);
%     a=load([GetComputerParams('BaseDir')   '/Code/SPM_precofigures/Segment.mat']);
%     a=load([getComputerParams('spmpreconfpath') 'Segment.mat']);
    a=load(['Segment.mat']);
    a.matlabbatch{1}.spm.spatial.preproc.data{1}=AFN; % Target
%     a.matlabbatch{1}.spm.spatial.preproc.opts.tpm=regexprep(a.matlabbatch{1}.spm.spatial.preproc.opts.tpm,'/local_data',GetComputerParams('BaseDir'));
%     a.matlabbatch{1}.spm.spatial.preproc.opts.tpm{1}='/data/spm8/tpm/grey.nii';
%     a.matlabbatch{1}.spm.spatial.preproc.opts.tpm{2}='/data/spm8/tpm/white.nii';
%     a.matlabbatch{1}.spm.spatial.preproc.opts.tpm{3}='/data/spm8/tpm/csf.nii';
    a.matlabbatch{1}.spm.spatial.preproc.opts.tpm=getComputerParams('tpm');
    if(nargin>4)
        if(~isempty(MaskFN))
            MFN=[TempPath  filesep 'MaskForSeg.nii'];
            copyfile(MaskFN,MFN,'f');
            a.matlabbatch{1}.spm.spatial.preproc.opts.msk{1}=MFN;
        end
    end
    spm_jobman('run',a.matlabbatch);
    TmpName= license('inuse');
    h=findobj('Name',['SPM8 (' TmpName(1).user '): Graphics']);
    close(h);
    mkdir(TrgP);
    delete([TrgP '*Seg*.*']);
    copyfile([TempPath  filesep '*.*'],TrgP,'f');
    % delete([GetComputerParams('BaseDir')   '/Temp/Coreg/*.*']);
    cd(CurP);
    rmdir(TempPath,'s');
    disp(['SPM Segmentation end for ' SFN ' ' datestr(now)]);
catch MExx
    cd(CurP);
    rmdir(TempPath,'s');
    error(['Error - SPM_Segment ' SFN]);
end