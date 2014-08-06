function TrgFN=n3_b1_clean_mask_FN(Src,FMaskFN,Force,TrgFN,FWHM,Dist)
% function TargetFN=n3_b1_clean(Src,DispProg)
if(~exist('Dist','var'))
    Dist=200;
end
if(isempty(Dist))
    Dist=200;
end
if(~exist('FWHM','var'))
    FWHM=0.15;
end
if(isempty(FWHM))
    FWHM=0.15;
end
if(~exist('Force','var'))
    Force=0;
end
if(isempty(Force))
    Force=0;
end
if(~exist('TrgFN','var'))
    TrgFN=[];
end
if(isempty(TrgFN))
    [A S]=fileparts(Src);
    TrgFN=[A filesep 'n3m_' S '.nii'];
end
if(exist(TrgFN,'file') && ~Force)
    disp(['Already exist - ' TrgFN]);
    return;
end

TmpName= license('inuse');
TmpName=TmpName(1).user;
TempP=[getComputerParams('temppath') '_N3_' TmpName '_' regexprep(datestr(now),'[-: ]','_') '_' num2str(floor(rand*1000)) filesep];
mkdir(TempP);

% TempP=[getComputerParams('temppath',1) filesep 'N3' filesep];
% delete([TempP '*.*']);
SFN=[TempP 'S.nii'];
MFN=[TempP 'M.nii'];
TFN=[TempP 'T.nii'];
SFNM=[TempP 'S.mnc'];
MFNM=[TempP 'M.mnc'];
TFNM=[TempP 'T.mnc'];
% copyfile(Src,SFN);
A=load_untouch_nii(Src);
A.img=A.img(end:-1:1,:,:);
save_untouch_nii(A,SFN);
% copyfile(Mask,MFN);
A=load_untouch_nii(FMaskFN);
A.img=A.img(end:-1:1,:,:);
save_untouch_nii(A,MFN);
% N3PathsStr='setenv PATH /u/peptibase3-ext/libermg1/OutsideCode/N3All/Out/Perl/site/bin:/u/peptibase3-ext/libermg1/OutsideCode/N3All/Out/Perl/bin:$PATH; setenv MANPATH /u/peptibase3-ext/libermg1/OutsideCode/N3All/Out/Perl/site/man:/u/peptibase3-ext/libermg1/OutsideCode/N3All/Out/Perl/man:$MANPATH; setenv PERL5LIB /u/peptibase3-ext/libermg1/OutsideCode/N3All/Out/MNIPerl/lib/; setenv PATH /u/peptibase3-ext/libermg1/OutsideCode/N3All/Out/bin:$PATH;/u/peptibase3-ext/libermg1/OutsideCode/N3All/Out/bin/';
N3PathsStr='';
system([N3PathsStr 'nii2mnc ' SFN]);
system([N3PathsStr 'nii2mnc ' MFN]);

disp('----1----');
% system(['nu_correct -clobber -shrink 2 -fwhm ' num2str(FWHM) ' -distance ' num2str(Dist) ' -mask ' MFNM ' ' SFNM ' ' TFNM]);
system([N3PathsStr 'nu_correct -clobber -shrink 2 -fwhm ' num2str(FWHM) ' -distance ' num2str(Dist) ' -mask ' MFNM ' ' SFNM ' ' TFNM]);
% system([N3PathsStr 'mnc2nii ' TFNM ' ' TFN]);
try
    [status, result] = system([N3PathsStr 'mnc2nii ' TFNM ' ' TempP 'a.nii']);
    B=load_untouch_nii([TempP 'a.nii']);
catch
    [status, result] = system([N3PathsStr 'mnc2nii ' TFNM ' ' TempP 'aa.nii']);
    B=load_untouch_nii([TempP 'aa.nii']);
end

A=load_untouch_nii(SFN);
% B=load_untouch_nii(TFN);
% delete([TempP 'a.nii']);
A.img=B.img;
save_untouch_nii(A,TrgFN);
rmdir(TempP,'s');