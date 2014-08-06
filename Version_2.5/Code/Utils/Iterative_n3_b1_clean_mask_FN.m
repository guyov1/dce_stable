function TrgFN=Iterative_n3_b1_clean_mask_FN(IterNum,Src,FMaskFN,Force,TrgFN,FWHM,Dist)
% function TrgFN=Iterative_n3_b1_clean_mask_FN(IterNum,Src,FMaskFN,Force,TrgFN,FWHM,Dist)
if(filesep=='\') % Windows
    TrgFN=Src;
    return;
end
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
if(~exist('TrgFN','var'))
    TrgFN=[];
end
if(isempty(TrgFN))
    [A S]=fileparts(Src);
    TrgFN=[A filesep 'n3m' num2str(IterNum) 'x_' S '.nii'];
end
if(~exist('Force','var'))
    Force=0;
end
if(isempty(Force))
    Force=0;
end
TargetFN=n3_b1_clean_mask_FN(Src,FMaskFN,Force,[],FWHM,Dist);
for i=2:IterNum-1
    TargetFN=n3_b1_clean_mask_FN(TargetFN,FMaskFN,Force,[],FWHM,Dist);
end
TrgFN=n3_b1_clean_mask_FN(TargetFN,FMaskFN,Force,TrgFN,FWHM,Dist);