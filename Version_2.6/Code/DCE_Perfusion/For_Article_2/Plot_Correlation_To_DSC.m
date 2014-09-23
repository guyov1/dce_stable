SliceNum      = 3;
Use_Mask      = true;
Use_Threshold = true;
Threshold_val = 200;

% Take values in the range 0.05 - 0.95 (to avoid margin effect)
MinRange = 0.20 ;
MaxRange = 0.80 ;

%% Create Normalized maps from DCE/DSC results
Path4Norm     = ['\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\' ...
                'Perfusion_DCE\coregFlow\DSC2DCE\'];
RefNii        = ['\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\' ...
                'Perfusion_DCE\coregFlow\DSC2DCE\Flow_Larsson_Relative_WM_30_6_brain_Thresholded_100_Normalized_0_1.nii'];
            FileNorm     = 'rdsc_oCBFlr';
NormalizeNii( Path4Norm, FileNorm, RefNii, Use_Threshold, Threshold_val)

%% Load data to correlate
Path_DCE      = ['\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\' ...
    'Perfusion_DCE\coregFlow\DSC2DCE\Flow_Larsson_Relative_WM_30_6_brain_Thresholded_200_Normalized_0_1.nii' ];
Path_DSC      = [ '\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\' ...
    'Perfusion_DCE\coregFlow\DSC2DCE\rdsc_oCBFlr_Thresholded_200_Normalized_0_1.nii'] ;

Path_Mask_1    = ['\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\forCorral\dce\ovl-wm-dce.nii'];
Path_Mask_2    = ['\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\forCorral\dce\ovl-wm-dce.nii'];
Path_Mask_3    = ['\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\forCorral\dce\ovl-wm-dce.nii'];
Path_Mask_4    = ['\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\forCorral\dce\ovl-bv-dce.nii'];

Mask1          = loadniidata(Path_Mask_1);
Mask2          = loadniidata(Path_Mask_2);
Mask3          = loadniidata(Path_Mask_3);
Mask4          = loadniidata(Path_Mask_4);
NiiDCE        = loadniidata(Path_DCE);
NiiDSC        = loadniidata(Path_DSC);

%% 

% figure;
% subplot(1,2,1);
% imshow(NiiDCE(:,:,SliceNum),'Colormap',jet(255));
% title(['DCE Map. Slice: ' num2str(SliceNum)]);
% subplot(1,2,2);
% imshow(NiiDSC(:,:,SliceNum),'Colormap',jet(255));
% title(['DSC Map. Slice: ' num2str(SliceNum)]);



CorrNiiMaps(Path_DCE, Path_DSC, SliceNum, MinRange, MaxRange, Mask1, Mask2, Mask3, Mask4, Use_Mask);