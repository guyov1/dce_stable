Name = 'rCBVCorrected';
Path = '\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\forCorral\dsc\DSCoMAN\';

A          = load_nii([Path Name '.nii']);
temp1      = load_nii([Path 'rCBVCorrected_flipped' '.nii']);

Nii2Flip   = A.img;
FlippedNii = flipdim(Nii2Flip,3);

%save_nii(B, [Path Name '_temp.nii']);

OutName    = [Path Name '_flipped.nii'];
refPath    = '\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\DCE-HTR\Results\2sec\Perfusion_DCE_with_delay';
refNiiFile = [refPath filesep 'Flow_Larsson.nii'];
Raw2Nii(FlippedNii, OutName, 'float32', refNiiFile);

temp1      = load_untouch_nii([Path 'DSC_4D' '.nii']);
temp2      = load_untouch_nii([Path 'rCBVCorrected_flipped' '.nii']);
