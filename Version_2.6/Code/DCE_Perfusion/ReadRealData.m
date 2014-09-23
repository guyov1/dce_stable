function [Subject_name, Subject_Path, WM_mask_absolute_path, Art_Mask, Vein_Mask, After_CTC_mat, DCECoregP] = ReadRealData()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Subject_name          = 'SmVl';
Subject_Path          = '\\fmri-t9\users\Moran\Stereotactic_Biopsy\PEKUROVSKI_NELENTINA\Study20140902_092528\DCE6min';
WM_mask_absolute_path = [Subject_Path  '\RefT1_WM_830.nii'];
%Art_Mask              = [Subject_Path  '\ManualArtMask.nii'];
Art_Mask              = [Subject_Path  '\InspectedRepVox.nii'];
Vein_Mask             = [Subject_Path  '\Veins_Mask.nii'];
After_CTC_mat         = [Subject_Path  '\AfterCTC.mat'];
%DCECoregP = [WorkingP 'DCEMainCoreged' filesep];
% \\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL
%DCECoregP             = [Subject_Path filesep 'DCE_out' filesep 'OrZe_20130811' filesep];
DCECoregP             = '\\fmri-t9\users\Moran\Stereotactic_Biopsy\PEKUROVSKI_NELENTINA\Study20140902_092528\DCE6min\PeNe_20140902_2sec\';



%Subject_name          = 'KUDISH_IRITH';

%Output_directory      = './Run_Output/';
%Subject_name          = 'SmVl';
% WM_mask_absolute_path = '\\fmri-t9\users\guyn\DCE_OUT\Guy_Test_1\SmVl_20120930\DCE4D_WM_Mask.nii';
% %After_CTC_mat         = '\\fmri-t9\users\Moran\DCE\TABASCO\KUDISH_IRITH\Study20131023_124549\DCE\dceout\KuIr_20131023\AfterCTC.mat';
% %AIFFindData_mat       = '\\fmri-t9\users\Moran\DCE\TABASCO\KUDISH_IRITH\Study20131023_124549\DCE\dceout\KuIr_20131023\AIFFindData.mat';
% After_CTC_mat         = 'D:\users\guyn\DCE_OUT\Guy_Test_1\SmVl_20120930\AfterCTC.mat';
% AIFFindData_mat       = 'D:\users\guyn\DCE_OUT\Guy_Test_1\SmVl_20120930\AIFFindData.mat';
% %Art_Mask              = '\\fmri-t9\users\Moran\DCE\TABASCO\KUDISH_IRITH\Study20131023_124549\DCE\ICAmasks\ARTcomponent.nii';
% %Vein_Mask             = '\\fmri-t9\users\Moran\DCE\TABASCO\KUDISH_IRITH\Study20131023_124549\DCE\ICAmasks\VEINScomponent.nii';
% Art_Mask              = '\\fmri-guy2\Dropbox\University\Msc\Thesis\SourceForge\DCE Data 2 sec\SMOLIAR_VLADIMIR\ICAmasks\ARTcomponent.nii';
% Vein_Mask             = '\\fmri-guy2\Dropbox\University\Msc\Thesis\SourceForge\DCE Data 2 sec\SMOLIAR_VLADIMIR\ICAmasks\VEINScomponent.nii';

% Subject_name          = 'MiAl';
% Subject_Path          = '\\fmri-t9\users\Moran\lesionVasClassification\GB\11RINDNER_GRETA_AUDREY\BL\';
% WM_mask_absolute_path = [Subject_Path 'prepfiles\wmref.nii'];
% Art_Mask              = [Subject_Path  'RiAu_20090518\manualArt.nii'];
% After_CTC_mat         = [Subject_Path  'RiAu_20090518\AfterCTC.mat'];

% Subject_name          = 'RoAs';
% Subject_Path          = '\\fmri-t9\users\Moran\lesionVasClassification\GB\01Rodity_Asaf\BL\';
% WM_mask_absolute_path = [Subject_Path 'prepfiles\wmref.nii'];
% Art_Mask              = [Subject_Path  'RoAs_20080122\manualArt.nii'];
% After_CTC_mat         = [Subject_Path  'RoAs_20080122\AfterCTC.mat'];

% Subject_name          = 'ZiYa';
% Subject_Path          = '\\FMRI-GUY2\Data\2 Sec\GlioBastoma\ZIFROT_YAAKOV\';
% WM_mask_absolute_path = [Subject_Path 'DCE\HTR\ZiYa_20130804\RefAuto_Base_WM_830.nii'];
% Art_Mask              = [Subject_Path  'ICAmasks\ARTcomponent.nii'];
% After_CTC_mat         = [Subject_Path  'DCE\HTR\ZiYa_20130804\AfterCTC.mat'];

% Subject_name          = 'ZiYa';
% Subject_Path          = '\\FMRI-GUY2\Data\2 Sec\GlioBastoma\ZIFROT_YAAKOV\';
% WM_mask_absolute_path = [Subject_Path 'DCE\STD\ZiYa_20130804\RefAuto_Base_WM_830.nii'];
% Art_Mask              = [Subject_Path  'ICAmasks\ARTcomponent.nii'];
% After_CTC_mat         = [Subject_Path  'DCE\STD\ZiYa_20130804\AfterCTC.mat'];

% Subject_name          = 'BaTa';
% Subject_Path          = '\\FMRI-GUY2\Data\2 Sec\Healthy + DSC\BARAK_TAL\';
% WM_mask_absolute_path = [Subject_Path 'DCE_out\BaTa_20131003x\RefAuto1_WM_830.nii'];
% Art_Mask              = [Subject_Path  'ICAmasks\ARTcomponent.nii'];
% After_CTC_mat         = [Subject_Path  'DCE_out\BaTa_20131003x\AfterCTC.mat'];

% Subject_name          = 'OrZe';
% Subject_Path          = '\\FMRI-GUY2\Data\2 Sec\Healthy\OR_DANA_ZEHAVA\';
% WM_mask_absolute_path = [Subject_Path 'DCE_out\OrZe_20130811\RefAuto1_WM_830.nii'];
% Art_Mask              = [Subject_Path  'ICAmasks\ARTcomponent.nii'];
% After_CTC_mat         = [Subject_Path  'DCE_out\OrZe_20130811\AfterCTC.mat'];

% Subject_name          = 'CoRa';
% Subject_Path          = '\\fmri-t9\users\Moran\Reports_for_deebi\019_COHEN_ZEDEK_RAHEL\Study20140402_100834_T1\';
% WM_mask_absolute_path = [Subject_Path  'DCE\CoRa_20140402\RefAuto1_WM_830.nii'];
% Art_Mask              = [Subject_Path  'DCE\ICA_masks\ARTcomponent.nii'];
% After_CTC_mat         = [Subject_Path  'DCE\CoRa_20140402\AfterCTC.mat'];

% Subject_name          = 'SmVl';
% Subject_Path          = '\\fmri-guy2\Data\2 Sec\GlioBastoma\SMOLIAR_VLADIMIR\';
% WM_mask_absolute_path = [Subject_Path 'DCE_GILAD_SmVl_20120930\RefAuto1_WM_830.nii'];
% Art_Mask              = [Subject_Path  'ICAmasks\ARTcomponent.nii'];
% After_CTC_mat         = [Subject_Path  'DCE_GILAD_SmVl_20120930\AfterCTC.mat'];

end

