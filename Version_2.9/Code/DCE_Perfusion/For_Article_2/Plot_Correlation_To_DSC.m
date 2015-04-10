close all;
SliceNum            = 4;
Use_Mask            = true;
Use_Threshold       = true;
Remove_Extreme_Vals = true;
Close_CorrPlot_Fig  = true;  % Don't display corr figures per segment
Remove_BV           = true;  % Ignore Blood Vessels Segment
Scatter_no_legend   = true; % Plot an additional scatter plot with no legend for article
Remove_Noisy        = false;  % Remove noisy pixels from correlation

font_size           = 35; % 25
font_legend         = 35; % 24
Threshold_val       = 350;

% Take values in the range 0.05 - 0.95 (to avoid margin effect)
MinRange = 0.00;
MaxRange = 1.00;
 
% Subj_Name             = 'REMEZ_YECHEZKEL';
% Subject_Path          = '\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\DCE-HTR';
% Path4Correl           = '\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\forCorral\dce';
%DCECoregP             = [Subject_Path  '\ReYe_20140615_2sec'];
Subj_Name             = 'PEKUROVSKI_NELENTINA';
Subject_Path          = '\\fmri-t9\users\Moran\DCE\HTR_GB\Stereotactic_Biopsy\PEKUROVSKI_NELENTINA\Study20140902_092528\DCE6min';
Path4Correl           = '\\fmri-t9\users\Moran\DCE\HTR_GB\Stereotactic_Biopsy\PEKUROVSKI_NELENTINA\forCorral\dce';
DCECoregP             = [Subject_Path  '\PeNe_20140902_2sec'];

WM_mask_absolute_path = [Subject_Path  '\RefT1_WM_830.nii'];
Art_Mask              = [Subject_Path  '\InspectedRepVox.nii'];
Vein_Mask             = [Subject_Path  '\Veins_Mask.nii'];
After_CTC_mat         = [Subject_Path  '\AfterCTC.mat'];
Brain_Extract_path    = [Subject_Path  '\Manual_BrainMask.nii'];
After_CTC_mat         = [Subject_Path  '\AfterCTC.mat'];

% Reference map to be able to write a new .nii output
RefNii                = [Path4Correl '\Flow_Larsson_Relative_WM_30_6.nii'];

% Maps name to check correlation
map_name        = 'Flow'; % Flow, MTT

if strcmp(map_name, 'Flow')
    FileNormDCE   = 'Flow_Larsson_Relative_WM_30_6';
    %FileNormDCE   = 'Flow_Larsson_model_Relative_WM_30_6';
    
    FileNormDSC   = 'rCBF';
elseif strcmp(map_name, 'MTT')
    FileNormDCE     = 'MTT_DCE_Sec_Brain_Extract';
    FileNormDSC     = 'rMTT';    
else
    error('-E- Unknown map requested!');
end

%% Remove noisy pixels
if Remove_Noisy
    
    if ~exist('CTC2D')
        display('-I- Loading CTC mat file...');
        load(After_CTC_mat);
    end
    CTC_4D             = Reshape2DCto4D(mat2cell(CTC2D,size(CTC2D,1),ones(1,size(CTC2D,2))),Msk2);
    mean_CTC           = mean(abs(CTC_4D),4);
    max_mean_CTC       = max(mean_CTC(:));
    med_mean_CTC       = median(mean_CTC(:));
    screen_val         = 0.001; % 0.00125, 0.001, 0.0015
    CTC_4D_Mask_By_Val = mean_CTC > screen_val*max_mean_CTC;
    
    % Brain Mask
    Brain_Mask_3D   = loadniidata(Brain_Extract_path);
    Brain_Mask_3D   = Brain_Mask_3D>0;
    
    num_before      = sum(sum(sum(Brain_Mask_3D>0)));
    display(['-I- Number of voxels before removing noisy voxels: ' num2str(num_before)]);
    
    % "And" the noisy mask with the brain mask
    Brain_Mask_3D   = Brain_Mask_3D .* CTC_4D_Mask_By_Val;
    
    num_after       = sum(sum(sum(Brain_Mask_3D>0)));
    display(['-I- Number of voxels after removing noisy voxels: ' num2str(num_after)]);
    
else
    % Brain Mask
    Brain_Mask_3D   = loadniidata(Brain_Extract_path);
    Brain_Mask_3D   = Brain_Mask_3D>0;
    
end

%% Create Normalized maps from DCE/DSC results

% Normalize DCE
Path_DCE     = NormalizeNii( Path4Correl, FileNormDCE, RefNii, Use_Threshold, Threshold_val);

% Normalize DSC
Path_DSC     = NormalizeNii( Path4Correl, FileNormDSC, RefNii, Use_Threshold, Threshold_val);


%% Load data to correlate
% Path_DCE      = ['\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\' ...
%     'Perfusion_DCE\coregFlow\DSC2DCE\Flow_Larsson_Relative_WM_30_6_brain_Thresholded_200_Normalized_0_1.nii' ];
% Path_DSC      = [ '\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\' ...
%     'Perfusion_DCE\coregFlow\DSC2DCE\rdsc_oCBFlr_Thresholded_200_Normalized_0_1.nii'] ;

Path_Mask_1    = [Path4Correl filesep 'ovl-wm-dce.nii'];
Path_Mask_2    = [Path4Correl filesep 'ovl-gm-dce.nii'];
Path_Mask_3    = [Path4Correl filesep 'ovl-lesion-dce.nii'];
Path_Mask_4    = [Path4Correl filesep 'ovl-bv-dce.nii'];
Mask_Names     = {{'WM'} {'GM'} {'Lesion'} {'bv'}};

Mask1          = loadniidata(Path_Mask_1);
Mask1          = Mask1 & Brain_Mask_3D;
Mask2          = loadniidata(Path_Mask_2);
Mask2          = Mask2 & Brain_Mask_3D;
Mask3          = loadniidata(Path_Mask_3);
Mask3          = Mask3 & Brain_Mask_3D;
Mask4          = loadniidata(Path_Mask_4);
Mask4          = Mask4 & Brain_Mask_3D;
NiiDCE         = loadniidata(Path_DCE);
NiiDSC         = loadniidata(Path_DSC);


%num_before      = sum(sum(sum(NiiDCE>0)));
%display(['-I- Number of voxels before removing noisy voxels from DCE: ' num2str(num_before)]);
NiiDCE(Brain_Mask_3D==0) = 0;
%num_after       = sum(sum(sum(NiiDCE>0)));
%display(['-I- Number of voxels after removing noisy voxels from DCE: ' num2str(num_after)]);
NiiDSC(Brain_Mask_3D==0) = 0;


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
            font_legend, map_name, Path4Correl, Close_CorrPlot_Fig, Remove_BV, Scatter_no_legend, Subj_Name);