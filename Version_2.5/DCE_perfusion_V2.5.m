%data location - Stable Version 2.1
% RUN - Test_On_Real_Data_Script

Subject_name          = 'CoRa';

Subject_Path          = '\\fmri-t9\users\Moran\Reports_for_deebi\019_COHEN_ZEDEK_RAHEL\Study20140402_100834_T1\';

WM_mask_absolute_path = [Subject_Path  'DCE\CoRa_20140402\RefAuto1_WM_830.nii'];

Art_Mask              = [Subject_Path  'DCE\ICA_masks\ARTcomponent.nii'];

After_CTC_mat         = [Subject_Path  'DCE\CoRa_20140402\AfterCTC.mat'];