function [DSCCoregP DSCCrgFNs] = DSCCoreg(DSCCoregP,src_dir)
% DSCCoregP=[dest_dir '\DSCMainCoreged\'];
mkdir(DSCCoregP);
DDSC=dir([src_dir filesep '*.nii']);
DSCFNs=strcat(src_dir,filesep,{DDSC.name})';
MatFN=RealignEstimate(DSCFNs);
DSCCrgFNs=CoregWrite(DSCFNs,MatFN,[],DSCCoregP);