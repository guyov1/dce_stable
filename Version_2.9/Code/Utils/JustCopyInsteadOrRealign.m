OrigP='\\fmri-t9\users\Moran\DCE\HTR_STROKE\REMEZ_YECHEZKEL\Study20140615_114415\DCE\ReYe_20140615\DCEMainNii\';
TrgP='\\fmri-t9\users\Moran\DCE\HTR_STROKE\REMEZ_YECHEZKEL\Study20140615_114415\DCE\ReYe_20140615\DCEMainCoreged\';
D=dir([OrigP '*.nii']);
for i=1:numel(D)
    copyfile([OrigP D(i).name],[TrgP 'Coreged_' D(i).name]);
end