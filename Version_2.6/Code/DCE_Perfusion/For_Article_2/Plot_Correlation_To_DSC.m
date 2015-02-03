close all;
SliceNum            = 4;
Use_Mask            = true;
Use_Threshold       = true;
Remove_Extreme_Vals = true;
Close_CorrPlot_Fig  = true; % Don't display corr figures per segment
Remove_BV           = true; % Ignore Blood Vessels Segment

font_size           = 25;
font_legend         = 24;
map_name            = 'Flow';

Threshold_val       = 350;

% Take values in the range 0.05 - 0.95 (to avoid margin effect)
MinRange = 0.00 ;
MaxRange = 1.00 ;

%% Create Normalized maps from DCE/DSC results

% Stroke Data
%'\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\forCorral\dce'
% Lesion Data
%'\\fmri-t9\users\Moran\Stereotactic_Biopsy\PEKUROVSKI_NELENTINA\forCorral\dce?ý'

Subj_Name     = 'REMEZ_YECHEZKEL';
Path4Norm     = ['\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\forCorral\dce'];
RefNii        = ['\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\forCorral\dce\Flow_Larsson_Relative_WM_30_6.nii'];

% Normalize DSC
%FileNorm     = 'rdsc_oCBFlr';
FileNorm     = 'rCBF';
%FileNorm     = 'rMTT';
Path_DSC = NormalizeNii( Path4Norm, FileNorm, RefNii, Use_Threshold, Threshold_val);
% Normalize DCE
%FileNorm     = 'DCE_MTT';
FileNorm     = 'Flow_Larsson_Relative_WM_30_6';
Path_DCE = NormalizeNii( Path4Norm, FileNorm, RefNii, Use_Threshold, Threshold_val);

%% Load data to correlate
% Path_DCE      = ['\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\' ...
%     'Perfusion_DCE\coregFlow\DSC2DCE\Flow_Larsson_Relative_WM_30_6_brain_Thresholded_200_Normalized_0_1.nii' ];
% Path_DSC      = [ '\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\' ...
%     'Perfusion_DCE\coregFlow\DSC2DCE\rdsc_oCBFlr_Thresholded_200_Normalized_0_1.nii'] ;

Path_Mask_1    = [Path4Norm filesep 'ovl-wm-dce.nii'];
Path_Mask_2    = [Path4Norm filesep 'ovl-gm-dce.nii'];
Path_Mask_3    = [Path4Norm filesep 'ovl-lesion-dce.nii'];
Path_Mask_4    = [Path4Norm filesep 'ovl-bv-dce.nii'];
Mask_Names     = {{'WM'} {'GM'} {'Lesion'} {'bv'}};

Mask1          = loadniidata(Path_Mask_1);
Mask2          = loadniidata(Path_Mask_2);
Mask3          = loadniidata(Path_Mask_3);
Mask4          = loadniidata(Path_Mask_4);
NiiDCE         = loadniidata(Path_DCE);
NiiDSC         = loadniidata(Path_DSC);

%% 

% figure;
% subplot(1,2,1);
% imshow(NiiDCE(:,:,SliceNum),'Colormap',jet(255));
% title(['DCE Map. Slice: ' num2str(SliceNum)]);
% subplot(1,2,2);
% imshow(NiiDSC(:,:,SliceNum),'Colormap',jet(255));
% title(['DSC Map. Slice: ' num2str(SliceNum)]);

CorrNiiMaps(Path_DCE, Path_DSC, SliceNum, MinRange, MaxRange, Mask1, Mask2, ...
            Mask3, Mask4, Use_Mask, Mask_Names, Remove_Extreme_Vals, font_size, ...
            font_legend, map_name, Path4Norm, Close_CorrPlot_Fig, Remove_BV, Subj_Name);