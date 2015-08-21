function [ ] = resultToNiiAndReport(resultStruct, time_vec_minutes, CTC2D, Chosen_AIF, Msk2, Brain_Mask_3D, Output_directory, DCECoregP, Sim_Struct, WM_mask_absolute_path)
%resultToNiiAndReport Summary of this function goes here




%% Take all results from Sim and return struct

useUptakeNoETM                         = Sim_Struct.useUptakeNoETM;
USE_ONE_GAUSSIAN                       = Sim_Struct.USE_ONE_GAUSSIAN;
USE_DOUBLE_GAUSSIAN                    = Sim_Struct.USE_DOUBLE_GAUSSIAN;
Ignore_Delay_Model_Selection           = Sim_Struct.Ignore_Delay_Model_Selection;
Data_Weight                            = Sim_Struct.Data_Weight;
AIC_Correction                         = Sim_Struct.AIC_Correction;
Threshold_Norm_Maps                    = Sim_Struct.Threshold_Norm_Maps;
Threshold_Val                          = Sim_Struct.Threshold_Val;

Flow_with_Delay                        = resultStruct.Flow_with_Delay;
Flow_no_Delay                          = resultStruct.Flow_no_Delay;
Delay_sec_by_Max_Val_with_Delay        = resultStruct.Delay_sec_by_Max_Val_with_Delay;
Delay_sec_by_Max_Val_no_Delay          = resultStruct.Delay_sec_by_Max_Val_no_Delay;
est_delay_by_AIF_correct               = resultStruct.est_delay_by_AIF_correct;
Est_IRF_with_Delay                     = resultStruct.Est_IRF_with_Delay;
Est_IRF_no_Delay                       = resultStruct.Est_IRF_no_Delay;
gaussian_param                         = resultStruct.gaussian_param;
fitted_gaussian                        = resultStruct.fitted_gaussian;
conv_result_Larss_with_Delay           = resultStruct.conv_result_Larss_with_Delay;
conv_result_Larss_no_Delay             = resultStruct.conv_result_Larss_no_Delay;
conv_result_Larss_no_Delay_High_F      = resultStruct.conv_result_Larss_no_Delay_High_F;
conv_result_Larss_no_Delay_no_E        = resultStruct.conv_result_Larss_no_Delay_no_E;
conv_result_Larss_no_Delay_no_E_High_F = resultStruct.conv_result_Larss_no_Delay_no_E_High_F;
conv_result_no_Delay_IRF               = resultStruct.conv_result_no_Delay_IRF;
conv_result_gaussian                   = resultStruct.conv_result_gaussian;
RMS_Larss_with_Delay                   = resultStruct.RMS_Larss_with_Delay;
RMS_Larss_no_Delay                     = resultStruct.RMS_Larss_no_Delay;
RMS_Larss_with_Delay_High_F            = resultStruct.RMS_Larss_with_Delay_High_F;
RMS_Larss_no_Delay_High_F              = resultStruct.RMS_Larss_no_Delay_High_F;
RMS_Larss_with_Delay_no_Ve             = resultStruct.RMS_Larss_with_Delay_no_Ve;
RMS_Larss_no_Delay_no_Ve               = resultStruct.RMS_Larss_no_Delay_no_Ve;
RMS_Larss_with_Delay_no_E              = resultStruct.RMS_Larss_with_Delay_no_E;
RMS_Larss_no_Delay_no_E                = resultStruct.RMS_Larss_no_Delay_no_E;
RMS_Larss_with_Delay_no_E_High_F       = resultStruct.RMS_Larss_with_Delay_no_E_High_F;
RMS_Larss_no_Delay_no_E_High_F         = resultStruct.RMS_Larss_no_Delay_no_E_High_F;
RMS_Larss_no_Delay_zero_params         = resultStruct.RMS_Larss_no_Delay_zero_params;
RMS_ht_no_Delay                        = resultStruct.RMS_ht_no_Delay;
RMS_gauss                              = resultStruct.RMS_gauss;
RMS_params_Gauss                       = resultStruct.RMS_params_Gauss;
fitted_double_gaussian                 = resultStruct.fitted_double_gaussian;
conv_result_double_gaussian            = resultStruct.conv_result_double_gaussian;
RMS_double_gauss                       = resultStruct.RMS_double_gauss;
RMS_params_double_gauss                = resultStruct.RMS_params_double_gauss;
Ktrans_with_Delay                      = resultStruct.Ktrans_with_Delay;
Ktrans_no_Delay                        = resultStruct.Ktrans_no_Delay;
E_with_Delay                           = resultStruct.E_with_Delay;
E_no_Delay                             = resultStruct.E_no_Delay;
E_with_Delay_no_Ve                     = resultStruct.E_with_Delay_no_Ve;
E_no_Delay_no_Ve                       = resultStruct.E_no_Delay_no_Ve;
Ktrans_with_Delay_High_F               = resultStruct.Ktrans_with_Delay_High_F;
Ktrans_no_Delay_High_F                 = resultStruct.Ktrans_no_Delay_High_F;
Vb_with_Delay                          = resultStruct.Vb_with_Delay;
Vb_no_Delay                            = resultStruct.Vb_no_Delay;
Vb_with_Delay_High_F                   = resultStruct.Vb_with_Delay_High_F;
Vb_no_Delay_High_F                     = resultStruct.Vb_no_Delay_High_F;
Vb_with_Delay_no_Ve                    = resultStruct.Vb_with_Delay_no_Ve;
Vb_no_Delay_no_Ve                      = resultStruct.Vb_no_Delay_no_Ve;
Vb_with_Delay_no_E                     = resultStruct.Vb_with_Delay_no_E;
Vb_no_Delay_no_E                       = resultStruct.Vb_no_Delay_no_E;
Vb_with_Delay_no_E_High_F              = resultStruct.Vb_with_Delay_no_E_High_F;
Vb_no_Delay_no_E_High_F                = resultStruct.Vb_no_Delay_no_E_High_F;
Ve_with_Delay                          = resultStruct.Ve_with_Delay;
Ve_no_Delay                            = resultStruct.Ve_no_Delay;
Ve_with_Delay_High_F                   = resultStruct.Ve_with_Delay_High_F;
Ve_no_Delay_High_F                     = resultStruct.Ve_no_Delay_High_F;
MTT_with_Delay                         = resultStruct.MTT_with_Delay;
MTT_no_Delay                           = resultStruct.MTT_no_Delay;
Ktrans_Patlak_with_Delay               = resultStruct.Ktrans_Patlak_with_Delay;
Ktrans_Patlak_no_Delay                 = resultStruct.Ktrans_Patlak_no_Delay;
Vb_Patlak_with_Delay                   = resultStruct.Vb_Patlak_with_Delay;
Vb_Patlak_no_Delay                     = resultStruct.Vb_Patlak_no_Delay;
MTT_Patlak_with_Delay                  = resultStruct.MTT_Patlak_with_Delay;
MTT_Patlak_no_Delay                    = resultStruct.MTT_Patlak_no_Delay;
double_gauss_param                    = resultStruct.double_gauss_param;

%% Reshape 2D to 3D and 4D

est_delay_by_AIF_correct_3D          = Reshape2DCto4D(est_delay_by_AIF_correct,Msk2);

t_delay_single_gauss_sec_3D          = Reshape2DCto4D(gaussian_param(1,:),Msk2);
sigma_seconds_single_gauss_3D        = Reshape2DCto4D(gaussian_param(2,:),Msk2);
Amp_single_gauss_3D                  = Reshape2DCto4D(gaussian_param(3,:),Msk2);

Est_IRF_with_Delay_4D                = Reshape2DCto4D(mat2cell(Est_IRF_with_Delay,size(Est_IRF_with_Delay,1),ones(1,size(Est_IRF_with_Delay,2))),Msk2);
Est_IRF_no_Delay_4D                  = Reshape2DCto4D(mat2cell(Est_IRF_no_Delay,size(Est_IRF_no_Delay,1),ones(1,size(Est_IRF_no_Delay,2))),Msk2);
CTC_4D                               = Reshape2DCto4D(mat2cell(CTC2D,size(CTC2D,1),ones(1,size(CTC2D,2))),Msk2);

conv_result_no_Delay_IRF_4D          = Reshape2DCto4D(mat2cell(conv_result_no_Delay_IRF,size(conv_result_no_Delay_IRF,1),ones(1,size(conv_result_no_Delay_IRF,2))),Msk2);
conv_result_gaussian_3D              = Reshape2DCto4D(mat2cell(conv_result_gaussian,size(conv_result_gaussian,1),ones(1,size(conv_result_gaussian,2))),Msk2);
RMS_ht_no_Delay_3D                   = Reshape2DCto4D(RMS_ht_no_Delay,Msk2);

RMS_Larss_with_Delay_3D              = Reshape2DCto4D(RMS_Larss_with_Delay,Msk2);
RMS_Larss_no_Delay_3D                = Reshape2DCto4D(RMS_Larss_no_Delay,Msk2);
RMS_Larss_with_Delay_High_F_3D       = Reshape2DCto4D(RMS_Larss_with_Delay_High_F,Msk2);
RMS_Larss_no_Delay_High_F_3D         = Reshape2DCto4D(RMS_Larss_no_Delay_High_F,Msk2);
RMS_Larss_with_Delay_no_Ve_3D        = Reshape2DCto4D(RMS_Larss_with_Delay_no_Ve,Msk2);
RMS_Larss_no_Delay_no_Ve_3D          = Reshape2DCto4D(RMS_Larss_no_Delay_no_Ve,Msk2);
RMS_Larss_with_Delay_no_E_3D         = Reshape2DCto4D(RMS_Larss_with_Delay_no_E,Msk2);
RMS_Larss_no_Delay_no_E_3D           = Reshape2DCto4D(RMS_Larss_no_Delay_no_E,Msk2);
RMS_Larss_with_Delay_no_E_High_F_3D  = Reshape2DCto4D(RMS_Larss_with_Delay_no_E_High_F,Msk2);
RMS_Larss_no_Delay_no_E_High_F_3D    = Reshape2DCto4D(RMS_Larss_no_Delay_no_E_High_F,Msk2);
RMS_Larss_no_Delay_zero_params_3D    = Reshape2DCto4D(RMS_Larss_no_Delay_zero_params,Msk2);

RMS_gauss_3D                         = Reshape2DCto4D(RMS_gauss,Msk2);
RMS_params_Gauss_3D                  = Reshape2DCto4D(RMS_params_Gauss,Msk2);
calculated_gaussian_3D               = Reshape2DCto4D(mat2cell(fitted_gaussian,size(fitted_gaussian,1),ones(1,size(fitted_gaussian,2))),Msk2);
calculated_double_gaussian_3D        = Reshape2DCto4D(mat2cell(fitted_double_gaussian,size(fitted_double_gaussian,1),ones(1,size(fitted_double_gaussian,2))),Msk2);
conv_result_double_gaussian_3D       = Reshape2DCto4D(mat2cell(conv_result_double_gaussian,size(conv_result_double_gaussian,1),ones(1,size(conv_result_double_gaussian,2))),Msk2);
% Larsson's
Flow_with_Delay_3D                   = Reshape2DCto4D(Flow_with_Delay,Msk2);
Flow_no_Delay_3D                     = Reshape2DCto4D(Flow_no_Delay,Msk2);
Delay_sec_by_Max_Val_with_Delay_3D   = Reshape2DCto4D(Delay_sec_by_Max_Val_with_Delay,Msk2);
Delay_sec_by_Max_Val_no_Delay_3D     = Reshape2DCto4D(Delay_sec_by_Max_Val_no_Delay,Msk2);
% Used AIF duplicated to all pixels
AIF_Larsson_Duplicated_3D            = repmat(Chosen_AIF,max(size(Est_IRF_no_Delay)),1);
AIF_Larsson_4D                       = Reshape2DCto4D(mat2cell(AIF_Larsson_Duplicated_3D,size(AIF_Larsson_Duplicated_3D,1),ones(1,size(AIF_Larsson_Duplicated_3D,2))),Msk2);

% Reshape double gaussian results
t_delay_1_double_gauss_seconds_3D    = Reshape2DCto4D(double_gauss_param(1,:),Msk2);
sigma_1_double_gauss_seconds_3D      = Reshape2DCto4D(double_gauss_param(2,:),Msk2);
amplitude_1_double_gauss_seconds_3D  = Reshape2DCto4D(double_gauss_param(3,:),Msk2);
t_delay_2_double_gauss_seconds_3D    = Reshape2DCto4D(double_gauss_param(4,:),Msk2);
sigma_2_double_gauss_seconds_3D      = Reshape2DCto4D(double_gauss_param(5,:),Msk2);
amplitude_2_double_gauss_seconds_3D  = Reshape2DCto4D(double_gauss_param(6,:),Msk2);
RMS_double_gauss_3D                  = Reshape2DCto4D(RMS_double_gauss,Msk2);
RMS_params_double_gauss_3D           = Reshape2DCto4D(RMS_params_double_gauss,Msk2);

Ktrans_with_Delay_3D                 = Reshape2DCto4D(Ktrans_with_Delay,Msk2);
Ktrans_no_Delay_3D                   = Reshape2DCto4D(Ktrans_no_Delay,Msk2);
Ktrans_with_Delay_High_F_3D          = Reshape2DCto4D(Ktrans_with_Delay_High_F,Msk2);
Ktrans_no_Delay_High_F_3D            = Reshape2DCto4D(Ktrans_no_Delay_High_F,Msk2);

E_with_Delay_3D                      = Reshape2DCto4D(E_with_Delay,Msk2);
E_no_Delay_3D                        = Reshape2DCto4D(E_no_Delay,Msk2);

E_with_Delay_no_Ve_3D                = Reshape2DCto4D(E_with_Delay_no_Ve,Msk2);
E_no_Delay_no_Ve_3D                  = Reshape2DCto4D(E_no_Delay_no_Ve,Msk2);

Vb_with_Delay_3D                     = Reshape2DCto4D(Vb_with_Delay,Msk2);
Vb_no_Delay_3D                       = Reshape2DCto4D(Vb_no_Delay,Msk2);
Vb_no_Delay_High_F_3D                = Reshape2DCto4D(Vb_no_Delay_High_F,Msk2);
Vb_with_Delay_High_F_3D              = Reshape2DCto4D(Vb_with_Delay_High_F,Msk2);

Vb_with_Delay_no_Ve_3D               = Reshape2DCto4D(Vb_with_Delay_no_Ve,Msk2);
Vb_no_Delay_no_Ve_3D                 = Reshape2DCto4D(Vb_no_Delay_no_Ve,Msk2);


Vb_no_Delay_no_E_3D                  = Reshape2DCto4D(Vb_no_Delay_no_E,Msk2);
Vb_with_Delay_no_E_3D                = Reshape2DCto4D(Vb_with_Delay_no_E,Msk2);
Vb_no_Delay_no_E_High_F_3D           = Reshape2DCto4D(Vb_no_Delay_no_E_High_F,Msk2);
Vb_with_Delay_no_E_High_F_3D         = Reshape2DCto4D(Vb_with_Delay_no_E_High_F,Msk2);

Ve_with_Delay_3D                     = Reshape2DCto4D(Ve_with_Delay,Msk2);
Ve_no_Delay_3D                       = Reshape2DCto4D(Ve_no_Delay,Msk2);
Ve_with_Delay_High_F_3D              = Reshape2DCto4D(Ve_with_Delay_High_F,Msk2);
Ve_no_Delay_High_F_3D                = Reshape2DCto4D(Ve_no_Delay_High_F,Msk2);

MTT_with_Delay_Min_3D                = Reshape2DCto4D(MTT_with_Delay,Msk2);
MTT_with_Delay_Sec_3D                = 60* MTT_with_Delay_Min_3D;
MTT_no_Delay_Min_3D                  = Reshape2DCto4D(MTT_no_Delay,Msk2);
MTT_no_Delay_Sec_3D                  = 60 * MTT_no_Delay_Min_3D;

MTT_with_Delay_Patlak_3D             = Reshape2DCto4D(MTT_Patlak_with_Delay,Msk2);
MTT_no_Delay_Patlak_3D               = Reshape2DCto4D(MTT_Patlak_no_Delay,Msk2);
Ktrans_Patlak_with_Delay_3D          = Reshape2DCto4D(Ktrans_Patlak_with_Delay,Msk2);
Ktrans_Patlak_no_Delay_3D            = Reshape2DCto4D(Ktrans_Patlak_no_Delay,Msk2);
Vb_Patlak_with_Delay_3D              = Reshape2DCto4D(Vb_Patlak_with_Delay,Msk2);
Vb_Patlak_no_Delay_3D                = Reshape2DCto4D(Vb_Patlak_no_Delay,Msk2);

%% PDF

% Title for PDF before displaying images
AddToLog(Output_directory,'idx_005','\\subsection*{Parameters Maps}');

% Reshape 2D to 4D
fig_num = figure;
subplot(1,2,1);
imagesc(mritransform(est_delay_by_AIF_correct_3D(:,:,2)));colorbar;%colormap(jet)
title('Slice #2 - Time Delay [Sec]');
subplot(1,2,2);
imagesc(mritransform(est_delay_by_AIF_correct_3D(:,:,3)));colorbar;
title('Slice #3 - Time Delay [Sec]');
gprint(fig_num,[Output_directory 'Time_Delay_Novel_Method.png']);
%gprint(fig_num,'Run_Output/Time_Delay_Gaussian.png');
AddToLog(Output_directory,'idx_006','TimeDelayNovel','Time_Delay_Novel_Method.png');

fig_num = figure;
subplot(1,2,1);
imagesc(mritransform(sigma_seconds_single_gauss_3D(:,:,2)));colorbar;
title('Slice #2 - Sigma [Sec]');
subplot(1,2,2);
imagesc(mritransform(sigma_seconds_single_gauss_3D(:,:,3)));colorbar;
title('Slice #3 - Sigma [Sec]');
gprint(fig_num,[Output_directory 'Sigma_Gaussian.png']);
%gprint(fig_num,'Run_Output/Sigma_Gaussian.png');
AddToLog(Output_directory,'idx_007','SigmaGaussian','Sigma_Gaussian.png');

fig_num = figure;
subplot(1,2,1);
imagesc(mritransform(RMS_params_Gauss_3D(:,:,2)));colorbar;
title('Slice #2 - RMS gaussian and ht');
subplot(1,2,2);
imagesc(mritransform(RMS_params_Gauss_3D(:,:,3)));colorbar;
title('Slice #3 - RMS gaussian and ht');
gprint(fig_num,[Output_directory 'RMS_gaussian_and_ht.png']);
%gprint(fig_num,'Run_Output/RMS_gaussian_and_ht.png');
AddToLog(Output_directory,'idx_008','RmsGaussianAndHt','RMS_gaussian_and_ht.png');

fig_num = figure;
subplot(1,2,1);
imagesc(mritransform(RMS_gauss_3D(:,:,2)));colorbar;
title('Slice #2 - RMS gauss to data');
subplot(1,2,2);
imagesc(mritransform(RMS_gauss_3D(:,:,3)));colorbar;
title('Slice #3 - RMS gauss to data');
gprint(fig_num,[Output_directory 'RMS_gauss_to_data.png']);
%gprint(fig_num,'Run_Output/RMS_gauss_to_data.png');
AddToLog(Output_directory,'idx_009','RmsGaussToData','RMS_gauss_to_data.png');

% Larsson's
fig_num = figure;
subplot(1,2,1);
imagesc(mritransform(Flow_no_Delay_3D(:,:,2)));colorbar;%colormap(jet)
title('Slice #2 - Flow - Larson [mL/100g/min]');
subplot(1,2,2);
imagesc(mritransform(Flow_no_Delay_3D(:,:,3)));colorbar;
title('Slice #3 - Flow - Larson [mL/100g/min]');
gprint(fig_num,[Output_directory 'Flow_Larsson.png']);
%gprint(fig_num,'Run_Output/Flow_Larsson.png');
AddToLog(Output_directory,'idx_010','FlowLarsson','Flow_Larsson.png');

fig_num = figure;
subplot(1,2,1);
imagesc(mritransform(Delay_sec_by_Max_Val_no_Delay_3D(:,:,2)));colorbar;%colormap(jet)
title('Slice #2 - Delay by Max Val- Larson [Sec]');
subplot(1,2,2);
imagesc(mritransform(Delay_sec_by_Max_Val_no_Delay_3D(:,:,3)));colorbar;
title('Slice #3 - Delay by Max Val - Larson [Sec]');
gprint(fig_num,[Output_directory 'Delay_Larsson.png']);
%gprint(fig_num,'Run_Output/Delay_Larsson.png');
AddToLog(Output_directory,'idx_011','DelayLarsson','Delay_Larsson.png');

DDCE      = dir([DCECoregP filesep '*.nii']);
% HDR File
DCEFNs    = strcat([DCECoregP filesep],{DDCE.name})';

%% Model Selection

% First maps should be Gilad's
nDataPoints      = length(time_vec_minutes);
zeros_map        = zeros(size(Flow_with_Delay_3D));
% Inf_map          = exp(100) * ones(size(Flow_Larsson_with_Delay_3D));

if Ignore_Delay_Model_Selection
    
    NParams          = [0 1 2 3 4];
    num_maps         = length(NParams);
    size_3D          = size(RMS_Larss_no_Delay_3D);
    RMSmaps          = zeros([size_3D num_maps]);
    
    RMSmaps(:,:,:,5) = RMS_Larss_with_Delay_3D;             % 5 Params
    %RMSmaps(:,:,:,8) = RMS_Larss_no_Delay_3D;               % 4 Params
    if useUptakeNoETM
        RMSmaps(:,:,:,4) = RMS_Larss_with_Delay_no_Ve_3D;      % 4 Params -----
    else
        RMSmaps(:,:,:,4) = RMS_Larss_with_Delay_High_F_3D;      % 4 Params -----
    end

    %RMSmaps(:,:,:,6) = RMS_Larss_no_Delay_High_F_3D;        % 3 Params
    RMSmaps(:,:,:,3) = RMS_Larss_with_Delay_no_E_3D;        % 3 Params -----
    %RMSmaps(:,:,:,4) = RMS_Larss_no_Delay_no_E_3D;          % 2 Params
    RMSmaps(:,:,:,2) = RMS_Larss_with_Delay_no_E_High_F_3D; % 2 Params -----
    %RMSmaps(:,:,:,2) = RMS_Larss_no_Delay_no_E_High_F_3D;   % 1 Params
    RMSmaps(:,:,:,1) = RMS_Larss_no_Delay_zero_params_3D;   % 0 Params
    
    [ ChosenByAIC_3D ]                      = Model_Selection( nDataPoints, RMSmaps, NParams, Data_Weight, AIC_Correction);
    
    % F
    F_Model_Selected_3D                     = zeros(size(Flow_with_Delay_3D));
    F_Model_Selected_3D (ChosenByAIC_3D==5) = Flow_with_Delay_3D     (ChosenByAIC_3D==5);
    if useUptakeNoETM
        F_Model_Selected_3D (ChosenByAIC_3D==4) = Flow_with_Delay_3D (ChosenByAIC_3D==4); % Put zero although its Inf
    else
        F_Model_Selected_3D (ChosenByAIC_3D==4) = zeros_map          (ChosenByAIC_3D==4); % Put zero although its Inf
    end
    F_Model_Selected_3D (ChosenByAIC_3D==3) = Flow_with_Delay_3D   (ChosenByAIC_3D==3);
    F_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2); % Put zero although its Inf
    F_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    % Vb
    Vb_Model_Selected_3D                    = zeros(size(Vb_no_Delay_3D));
    Vb_Model_Selected_3D(ChosenByAIC_3D==5) = Vb_with_Delay_3D             (ChosenByAIC_3D==5);
    if useUptakeNoETM
        Vb_Model_Selected_3D(ChosenByAIC_3D==4) = Vb_with_Delay_no_Ve_3D      (ChosenByAIC_3D==4);
    else
        Vb_Model_Selected_3D(ChosenByAIC_3D==4) = Vb_with_Delay_High_F_3D      (ChosenByAIC_3D==4);
    end
    Vb_Model_Selected_3D(ChosenByAIC_3D==3) = Vb_with_Delay_no_E_3D        (ChosenByAIC_3D==3);
    Vb_Model_Selected_3D(ChosenByAIC_3D==2) = Vb_with_Delay_no_E_High_F_3D (ChosenByAIC_3D==2);
    Vb_Model_Selected_3D(ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    % Ve
    Ve_Model_Selected_3D                    = zeros(size(Ve_no_Delay_3D));
    Ve_Model_Selected_3D(ChosenByAIC_3D==5) = Ve_with_Delay_3D             (ChosenByAIC_3D==5);
    if useUptakeNoETM
        Ve_Model_Selected_3D(ChosenByAIC_3D==4) = zeros_map                (ChosenByAIC_3D==4);
    else
        Ve_Model_Selected_3D(ChosenByAIC_3D==4) = Ve_with_Delay_High_F_3D  (ChosenByAIC_3D==4);
    end
    
    Ve_Model_Selected_3D(ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
    Ve_Model_Selected_3D(ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
    Ve_Model_Selected_3D(ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    % E
    E_Model_Selected_3D                     = zeros(size(E_no_Delay_3D));
    E_Model_Selected_3D (ChosenByAIC_3D==5) = E_with_Delay_3D              (ChosenByAIC_3D==5);
    if useUptakeNoETM
        E_Model_Selected_3D (ChosenByAIC_3D==4) = E_with_Delay_no_Ve_3D    (ChosenByAIC_3D==4); % Ktrans is not E
    else
        E_Model_Selected_3D (ChosenByAIC_3D==4) = zeros_map                    (ChosenByAIC_3D==4); % Ktrans is not E
    end
    E_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
    E_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
    E_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    % Ktrans
    Ktrans_Model_Selected_3D                     = zeros(size(Ktrans_no_Delay_3D));
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==5) = Ktrans_with_Delay_3D         (ChosenByAIC_3D==5);
    if useUptakeNoETM
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==4) = E_with_Delay_no_Ve_3D  (ChosenByAIC_3D==4) .* Flow_with_Delay_3D  (ChosenByAIC_3D==4);
    else
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==4) = Ktrans_with_Delay_High_F_3D  (ChosenByAIC_3D==4);
    end
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    % Model Selection Map
    MeanFN=[Output_directory 'Model_Selection_Map.nii'];
    Raw2Nii(ChosenByAIC_3D,MeanFN,'float32',DCEFNs{1});
    
else % Include Delay
    %%
    NParams          = [0 1 2 2 3 3 4 4 5];
    num_maps         = length(NParams);
    size_3D          = size(RMS_Larss_no_Delay_3D);
    RMSmaps          = zeros([size_3D num_maps]);
    
    RMSmaps(:,:,:,9) = RMS_Larss_with_Delay_3D;             % 5 Params
    RMSmaps(:,:,:,8) = RMS_Larss_no_Delay_3D;               % 4 Params
    if useUptakeNoETM
        RMSmaps(:,:,:,7) = RMS_Larss_with_Delay_no_Ve_3D;      % 4 Params -----
        RMSmaps(:,:,:,6) = RMS_Larss_no_Delay_no_Ve_3D;        % 3 Params
    else
        RMSmaps(:,:,:,7) = RMS_Larss_with_Delay_High_F_3D;      % 4 Params -----
        RMSmaps(:,:,:,6) = RMS_Larss_no_Delay_High_F_3D;        % 3 Params
    end
    RMSmaps(:,:,:,5) = RMS_Larss_with_Delay_no_E_3D;        % 3 Params -----
    RMSmaps(:,:,:,4) = RMS_Larss_no_Delay_no_E_3D;          % 2 Params
    RMSmaps(:,:,:,3) = RMS_Larss_with_Delay_no_E_High_F_3D; % 2 Params -----
    RMSmaps(:,:,:,2) = RMS_Larss_no_Delay_no_E_High_F_3D;   % 1 Params
    RMSmaps(:,:,:,1) = RMS_Larss_no_Delay_zero_params_3D;   % 0 Params
    
    [ ChosenByAIC_3D ]                      = Model_Selection( nDataPoints, RMSmaps, NParams, Data_Weight, AIC_Correction);
    
    % Delay
    AIF_Delay_Model_Selected_3D                    = zeros(size(est_delay_by_AIF_correct_3D));
    AIF_Delay_Model_Selected_3D(ChosenByAIC_3D==9) = est_delay_by_AIF_correct_3D (ChosenByAIC_3D==9);
    AIF_Delay_Model_Selected_3D(ChosenByAIC_3D==7) = est_delay_by_AIF_correct_3D (ChosenByAIC_3D==7);
    AIF_Delay_Model_Selected_3D(ChosenByAIC_3D==5) = est_delay_by_AIF_correct_3D (ChosenByAIC_3D==5);
    AIF_Delay_Model_Selected_3D(ChosenByAIC_3D==3) = est_delay_by_AIF_correct_3D (ChosenByAIC_3D==3);
    
    % F
    F_Model_Selected_3D                     = zeros(size(Flow_with_Delay_3D));
    F_Model_Selected_3D (ChosenByAIC_3D==9) = Flow_with_Delay_3D   (ChosenByAIC_3D==9);
    F_Model_Selected_3D (ChosenByAIC_3D==8) = Flow_no_Delay_3D     (ChosenByAIC_3D==8);
    if useUptakeNoETM
        F_Model_Selected_3D (ChosenByAIC_3D==7) = Flow_with_Delay_3D           (ChosenByAIC_3D==7); % Put zero although its Inf
        F_Model_Selected_3D (ChosenByAIC_3D==6) = Flow_no_Delay_3D             (ChosenByAIC_3D==6); % Put zero although its Inf
    else
        F_Model_Selected_3D (ChosenByAIC_3D==7) = zeros_map                    (ChosenByAIC_3D==7); % Put zero although its Inf
        F_Model_Selected_3D (ChosenByAIC_3D==6) = zeros_map                    (ChosenByAIC_3D==6); % Put zero although its Inf
    end
    F_Model_Selected_3D (ChosenByAIC_3D==5) = Flow_with_Delay_3D   (ChosenByAIC_3D==5);
    F_Model_Selected_3D (ChosenByAIC_3D==4) = Flow_no_Delay_3D     (ChosenByAIC_3D==4);
    F_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3); % Put zero although its Inf
    F_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2); % Put zero although its Inf
    F_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    
    % Vb
    Vb_Model_Selected_3D                    = zeros(size(Vb_no_Delay_3D));
    Vb_Model_Selected_3D(ChosenByAIC_3D==9) = Vb_with_Delay_3D             (ChosenByAIC_3D==9);
    Vb_Model_Selected_3D(ChosenByAIC_3D==8) = Vb_no_Delay_3D               (ChosenByAIC_3D==8);
    if useUptakeNoETM
        Vb_Model_Selected_3D(ChosenByAIC_3D==7) = Vb_with_Delay_no_Ve_3D   (ChosenByAIC_3D==7);
        Vb_Model_Selected_3D(ChosenByAIC_3D==6) = Vb_no_Delay_no_Ve_3D     (ChosenByAIC_3D==6);
    else
        Vb_Model_Selected_3D(ChosenByAIC_3D==7) = Vb_with_Delay_High_F_3D  (ChosenByAIC_3D==7);
        Vb_Model_Selected_3D(ChosenByAIC_3D==6) = Vb_no_Delay_High_F_3D    (ChosenByAIC_3D==6);
    end
    Vb_Model_Selected_3D(ChosenByAIC_3D==5) = Vb_with_Delay_no_E_3D        (ChosenByAIC_3D==5);
    Vb_Model_Selected_3D(ChosenByAIC_3D==4) = Vb_no_Delay_no_E_3D          (ChosenByAIC_3D==4);
    Vb_Model_Selected_3D(ChosenByAIC_3D==3) = Vb_with_Delay_no_E_High_F_3D (ChosenByAIC_3D==3);
    Vb_Model_Selected_3D(ChosenByAIC_3D==2) = Vb_no_Delay_no_E_High_F_3D   (ChosenByAIC_3D==2);
    Vb_Model_Selected_3D(ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    % Ve
    Ve_Model_Selected_3D                    = zeros(size(Ve_no_Delay_3D));
    Ve_Model_Selected_3D(ChosenByAIC_3D==9) = Ve_with_Delay_3D             (ChosenByAIC_3D==9);
    Ve_Model_Selected_3D(ChosenByAIC_3D==8) = Ve_no_Delay_3D               (ChosenByAIC_3D==8);
    if useUptakeNoETM
        Ve_Model_Selected_3D(ChosenByAIC_3D==7) = zeros_map      (ChosenByAIC_3D==7);
        Ve_Model_Selected_3D(ChosenByAIC_3D==6) = zeros_map        (ChosenByAIC_3D==6);
    else
        Ve_Model_Selected_3D(ChosenByAIC_3D==7) = Ve_with_Delay_High_F_3D      (ChosenByAIC_3D==7);
        Ve_Model_Selected_3D(ChosenByAIC_3D==6) = Ve_no_Delay_High_F_3D        (ChosenByAIC_3D==6);
    end
    Ve_Model_Selected_3D(ChosenByAIC_3D==5) = zeros_map                    (ChosenByAIC_3D==5);
    Ve_Model_Selected_3D(ChosenByAIC_3D==4) = zeros_map                    (ChosenByAIC_3D==4);
    Ve_Model_Selected_3D(ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
    Ve_Model_Selected_3D(ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
    Ve_Model_Selected_3D(ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    % E
    E_Model_Selected_3D                     = zeros(size(E_no_Delay_3D));
    E_Model_Selected_3D (ChosenByAIC_3D==9) = E_with_Delay_3D              (ChosenByAIC_3D==9);
    E_Model_Selected_3D (ChosenByAIC_3D==8) = E_no_Delay_3D                (ChosenByAIC_3D==8);
    if useUptakeNoETM
        E_Model_Selected_3D (ChosenByAIC_3D==7) = E_with_Delay_no_Ve_3D    (ChosenByAIC_3D==7); % Ktrans is not E
        E_Model_Selected_3D (ChosenByAIC_3D==6) = E_no_Delay_no_Ve_3D      (ChosenByAIC_3D==6); % Ktrans is not E
    else
        E_Model_Selected_3D (ChosenByAIC_3D==7) = zeros_map                    (ChosenByAIC_3D==7); % Ktrans is not E
        E_Model_Selected_3D (ChosenByAIC_3D==6) = zeros_map                    (ChosenByAIC_3D==6); % Ktrans is not E
    end
    E_Model_Selected_3D (ChosenByAIC_3D==5) = zeros_map                    (ChosenByAIC_3D==5);
    E_Model_Selected_3D (ChosenByAIC_3D==4) = zeros_map                    (ChosenByAIC_3D==4);
    E_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
    E_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
    E_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
    
    % Ktrans
    Ktrans_Model_Selected_3D                     = zeros(size(Ktrans_no_Delay_3D));
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==9) = Ktrans_with_Delay_3D         (ChosenByAIC_3D==9);
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==8) = Ktrans_no_Delay_3D           (ChosenByAIC_3D==8);
    if useUptakeNoETM
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==7) = E_with_Delay_no_Ve_3D  (ChosenByAIC_3D==7) .*  Flow_with_Delay_3D  (ChosenByAIC_3D==7);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==6) = E_no_Delay_no_Ve_3D    (ChosenByAIC_3D==6) .*  Flow_no_Delay_3D    (ChosenByAIC_3D==6);
    else
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==7) = Ktrans_with_Delay_High_F_3D  (ChosenByAIC_3D==7);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==6) = Ktrans_no_Delay_High_F_3D    (ChosenByAIC_3D==6);
    end
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==5) = zeros_map                    (ChosenByAIC_3D==5);
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==4) = zeros_map                    (ChosenByAIC_3D==4);
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
    Ktrans_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);

    % Model Selection Map
    MeanFN=[Output_directory 'Model_Selection_Map.nii'];
    Raw2Nii(ChosenByAIC_3D,MeanFN,'float32',DCEFNs{1});
    %%
end



%% Maps

if (USE_ONE_GAUSSIAN)
    MeanFN=[Output_directory 'Delay_Single_Gaussian.nii'];
    Raw2Nii(t_delay_single_gauss_sec_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'Sigma_Single_Gaussian.nii'];
    Raw2Nii(sigma_seconds_single_gauss_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'Single_Gaussian_Amplitude.nii'];
    Raw2Nii(Amp_single_gauss_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'Est_Single_Gaussian_Filter.nii'];
    Raw2Nii(calculated_gaussian_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'Conv_Res_with_Single_Gaussian.nii'];
    Raw2Nii(conv_result_gaussian_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'RMS_Res_with_Single_Gauss_and_Ct.nii'];
    Raw2Nii(RMS_gauss_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'RMS_Ht_with_Fitted_Single_Gaussian.nii'];
    Raw2Nii(RMS_params_Gauss_3D,MeanFN,'float32',DCEFNs{1});
    
end

if (USE_DOUBLE_GAUSSIAN)
    
    MeanFN=[Output_directory 'Time_Delay_1_double_gaussian.nii'];
    Raw2Nii(t_delay_1_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'Sigma_1_double_gaussian.nii'];
    Raw2Nii(sigma_1_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'Amplitude_1_double_gaussian.nii'];
    Raw2Nii(amplitude_1_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'Time_Delay_2_double_gaussian.nii'];
    Raw2Nii(t_delay_2_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'Sigma_2_double_gaussian.nii'];
    Raw2Nii(sigma_2_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'Amplitude_2_double_gaussian.nii'];
    Raw2Nii(amplitude_2_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'Est_Double_Gaussian_Filter.nii'];
    Raw2Nii(calculated_double_gaussian_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'Conv_Res_with_Double_Gaussian.nii'];
    Raw2Nii(conv_result_double_gaussian_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'RMS_Res_with_Double_Gauss_and_Ct.nii'];
    Raw2Nii(RMS_double_gauss_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[Output_directory 'RMS_Ht_with_Fitted_Double_Gaussian.nii'];
    Raw2Nii(RMS_params_double_gauss_3D,MeanFN,'float32',DCEFNs{1});
    
end


%% Create a mask according to CTC median value

% Get the median of the entire volumes
% mean_CTC                             = mean(abs(CTC_4D),4);
% max_mean_CTC                         = max(mean_CTC(:));
% med_mean_CTC                         = median(mean_CTC(:));
median_CTC                           = median(abs(CTC_4D), 4);
max_median_CTC                       = max(median_CTC(:));
screen_val                           = 0.15 / 100; % Percent threshold

%CTC_4D_Mask_By_Val                   = mean_CTC > screen_val*max_mean_CTC;
CTC_4D_Mask_By_Val                   = median_CTC >= screen_val*max_median_CTC;
CTC_4D_Mask_negative_By_Val          = median_CTC <  screen_val*max_median_CTC;

% ------------------- Delay -----------------------------

%% Screen noisy voxels
% Get the median of the entire volumes
% MdA           = median(abs(CTC_4D), 4);
% The BAT map to be masked
MaskedBAT     = est_delay_by_AIF_correct_3D;
% Threshold for median value
% MdA_threshold = 0.00007;
% % Zero all valus below the threshold
% MaskedBAT(MdA < MdA_threshold) = 0;

% Zero all valus below the threshold
MaskedBAT(CTC_4D_Mask_negative_By_Val) = 0;

% Gilad's debug test to see the resulted image
tmp = squeeze(est_delay_by_AIF_correct_3D(:,:,4)); tmp(CTC_4D_Mask_negative_By_Val(:,:,4))= NaN;figure;imagesc(rot90(tmp))

% Output the map
Raw2Nii(MaskedBAT, [Output_directory 'BAT_ACOPED_Noise_Masked.nii'],'float32',DCEFNs{1});

% Same with model selection
MaskedBAT_model_selection                    = zeros(size(est_delay_by_AIF_correct_3D));
MaskedBAT_model_selection(ChosenByAIC_3D==9) = MaskedBAT (ChosenByAIC_3D==9);
MaskedBAT_model_selection(ChosenByAIC_3D==7) = MaskedBAT (ChosenByAIC_3D==7);
MaskedBAT_model_selection(ChosenByAIC_3D==5) = MaskedBAT (ChosenByAIC_3D==5);
MaskedBAT_model_selection(ChosenByAIC_3D==3) = MaskedBAT (ChosenByAIC_3D==3);

Raw2Nii(MaskedBAT_model_selection,[Output_directory 'BAT_ACOPED_Noise_Masked_Model_Selection.nii'],'float32',DCEFNs{1});

%%

MeanFN=[Output_directory 'BAT_ACOPED.nii'];
Raw2Nii(est_delay_by_AIF_correct_3D,MeanFN,'float32',DCEFNs{1});

if ~Ignore_Delay_Model_Selection
    MeanFN=[Output_directory 'BAT_ACOPED_Model_Selection.nii'];
    %%
    MeanFN=[Output_directory 'BAT_ACOPED_Model_Selection_Test1.nii'];
    Raw2Nii(AIF_Delay_Model_Selected_3D,MeanFN,'float32',DCEFNs{1});
    %%
end

MeanFN=[Output_directory 'Delay_by_Max_Val_Larsson_no_Delay.nii'];
Raw2Nii(Delay_sec_by_Max_Val_no_Delay_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Delay_by_Max_Val_Larsson_with_Delay.nii'];
Raw2Nii(Delay_sec_by_Max_Val_with_Delay_3D,MeanFN,'float32',DCEFNs{1});

% -------------------- F ----------------------------

MeanFN=[Output_directory 'Flow_Larsson_with_Delay.nii'];
Raw2Nii(Flow_with_Delay_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Flow_Larsson_no_Delay.nii'];
Raw2Nii(Flow_no_Delay_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Flow_Model_Selection.nii'];
Raw2Nii(F_Model_Selected_3D,MeanFN,'float32',DCEFNs{1});

% -------------------- E ----------------------------

MeanFN=[Output_directory 'E_with_Delay.nii'];
Raw2Nii(E_with_Delay_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'E_no_Delay.nii'];
Raw2Nii(E_no_Delay_3D,MeanFN,'float32',DCEFNs{1});

%MeanFN=[Output_directory 'E_with_Delay_High_F.nii'];
%Raw2Nii(Ktrans_with_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});

%MeanFN=[Output_directory 'E_no_Delay_High_F.nii'];
%Raw2Nii(Ktrans_no_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'E_Model_Selected.nii'];
Raw2Nii(E_Model_Selected_3D,MeanFN,'float32',DCEFNs{1});

% -------------------- Ktrans ----------------------------

MeanFN=[Output_directory 'Ktrans_with_Delay.nii'];
Raw2Nii(Ktrans_with_Delay_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Ktrans_no_Delay.nii'];
Raw2Nii(Ktrans_no_Delay_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Ktrans_with_Delay_High_F.nii'];
Raw2Nii(Ktrans_with_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Ktrans_no_Delay_High_F.nii'];
Raw2Nii(Ktrans_no_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Ktrans_with_Delay_Patlak.nii'];
Raw2Nii(Ktrans_Patlak_with_Delay_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Ktrans_no_Delay_Patlak.nii'];
Raw2Nii(Ktrans_Patlak_no_Delay_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Ktrans_Model_Selected.nii'];
Raw2Nii(Ktrans_Model_Selected_3D,MeanFN,'float32',DCEFNs{1});

% -------------------- Vb ----------------------------

MeanFN=[Output_directory 'Vb_with_Delay.nii'];
Raw2Nii(Vb_with_Delay_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Vb_no_Delay.nii'];
Raw2Nii(Vb_no_Delay_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Vb_with_Delay_High_F.nii'];
Raw2Nii(Vb_with_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Vb_no_Delay_High_F.nii'];
Raw2Nii(Vb_no_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Vb_with_Delay_no_E.nii'];
Raw2Nii(Vb_with_Delay_no_E_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Vb_no_Delay_no_E.nii'];
Raw2Nii(Vb_no_Delay_no_E_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Vb_with_Delay_no_E_High_F.nii'];
Raw2Nii(Vb_with_Delay_no_E_High_F_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Vb_no_Delay_no_E_High_F.nii'];
Raw2Nii(Vb_no_Delay_no_E_High_F_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Vb_Model_Selection.nii'];
Raw2Nii(Vb_Model_Selected_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Vb_with_Delay_Patlak.nii'];
Raw2Nii(Vb_Patlak_with_Delay_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Vb_no_Delay_Patlak.nii'];
Raw2Nii(Vb_Patlak_no_Delay_3D,MeanFN,'float32',DCEFNs{1});

% -------------------- Ve ----------------------------

MeanFN=[Output_directory 'Ve_with_Delay.nii'];
Raw2Nii(Ve_with_Delay_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Ve_no_Delay.nii'];
Raw2Nii(Ve_no_Delay_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Ve_with_Delay_High_F.nii'];
Raw2Nii(Ve_with_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Ve_no_Delay_High_F.nii'];
Raw2Nii(Ve_no_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Ve_Model_Selection.nii'];
Raw2Nii(Ve_Model_Selected_3D,MeanFN,'float32',DCEFNs{1});

% -------------------- MTT ----------------------------

MeanFN=[Output_directory 'MTT_with_Delay_Sec.nii'];
Raw2Nii(MTT_with_Delay_Sec_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'MTT_no_Delay_Sec.nii'];
Raw2Nii(MTT_no_Delay_Sec_3D,MeanFN,'float32',DCEFNs{1});

MTT_no_Delay_Sec_3D_brain_extract                    = zeros(size(MTT_no_Delay_Sec_3D));
MTT_no_Delay_Sec_3D_brain_extract(Brain_Mask_3D > 0) = MTT_no_Delay_Sec_3D(Brain_Mask_3D > 0);
MeanFN=[Output_directory 'MTT_no_Delay_Sec_Brain_Extract.nii'];
Raw2Nii(MTT_no_Delay_Sec_3D_brain_extract,MeanFN,'float32',DCEFNs{1});


MTT_with_Delay_Sec_3D_Norm_0_1                                     = MTT_with_Delay_Sec_3D ./ max(max(max(MTT_with_Delay_Sec_3D)));
MeanFN=[Output_directory 'MTT_no_Delay_Sec_Norm_0_1.nii'];
Raw2Nii(MTT_with_Delay_Sec_3D_Norm_0_1,MeanFN,'float32',DCEFNs{1});


max_idx = length(time_vec_minutes); %130
min_idx = 1;
est_MTT_noise = cumtrapz(time_vec_minutes(min_idx:max_idx),Est_IRF_with_Delay_4D(:,:,:,min_idx:max_idx),4); est_MTT_noise = est_MTT_noise(:,:,:,end);

MeanFN=[Output_directory 'MTT_Test.nii'];
Raw2Nii(est_MTT_noise*60,MeanFN,'float32',DCEFNs{1});


MeanFN=[Output_directory 'MTT_with_Delay_Patlak.nii'];
Raw2Nii(MTT_with_Delay_Patlak_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'MTT_no_Delay_Patlak.nii'];
Raw2Nii(MTT_no_Delay_Patlak_3D,MeanFN,'float32',DCEFNs{1});

% -------------------- RMS ----------------------------

MeanFN=[Output_directory 'Log_RMS_Conv_Res_with_IRF_and_Ct.nii'];
Raw2Nii(log(RMS_ht_no_Delay_3D),MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Log_RMS_Larss_with_Delay.nii'];
Raw2Nii(log(RMS_Larss_with_Delay_3D),MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Log_RMS_Larss_No_Delay.nii'];
Raw2Nii(log(RMS_Larss_no_Delay_3D),MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Log_RMS_Larss_with_Delay_High_F.nii'];
Raw2Nii(log(RMS_Larss_with_Delay_High_F_3D),MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'RMS_Larss_with_Delay_no_Ve_3D.nii'];
Raw2Nii(log(RMS_Larss_with_Delay_no_Ve_3D),MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'RMS_Larss_no_Delay_no_Ve_3D.nii'];
Raw2Nii(log(RMS_Larss_no_Delay_no_Ve_3D),MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Log_RMS_Larss_no_Delay_High_F.nii'];
Raw2Nii(log(RMS_Larss_no_Delay_High_F_3D),MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Log_RMS_Larss_with_Delay_no_E.nii'];
Raw2Nii(log(RMS_Larss_with_Delay_no_E_3D),MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Log_RMS_Larss_no_Delay_no_E.nii'];
Raw2Nii(log(RMS_Larss_no_Delay_no_E_3D),MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Log_RMS_Larss_with_Delay_no_E_High_F.nii'];
Raw2Nii(log(RMS_Larss_with_Delay_no_E_High_F_3D),MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Log_RMS_Larss_no_Delay_no_E_High_F.nii'];
Raw2Nii(log(RMS_Larss_no_Delay_no_E_High_F_3D),MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Log_RMS_Larss_no_Delay_zero_params.nii'];
Raw2Nii(log(RMS_Larss_no_Delay_zero_params_3D),MeanFN,'float32',DCEFNs{1});

% -------------------- General ----------------------------

MeanFN=[Output_directory 'AIF_Used_Larsson.nii'];
Raw2Nii(AIF_Larsson_4D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'CTC_per_voxel.nii'];
Raw2Nii(CTC_4D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'CTC_Mask_per_voxel_' num2str(screen_val*100) '.nii'];
Raw2Nii(CTC_4D_Mask_By_Val,MeanFN,'float32',DCEFNs{1});

CTC_4D_Mask_By_Val_brain_exract                         = zeros(size(CTC_4D_Mask_By_Val));
CTC_4D_Mask_By_Val_brain_exract(Brain_Mask_3D > 0) = CTC_4D_Mask_By_Val(Brain_Mask_3D > 0);

MeanFN=[Output_directory 'CTC_Mask_per_voxel_' num2str(screen_val*100) '_Brain_Extract.nii'];
Raw2Nii(CTC_4D_Mask_By_Val_brain_exract,MeanFN,'float32',DCEFNs{1});


est_delay_by_AIF_correct_no_noisy_3D = zeros(size(est_delay_by_AIF_correct_3D));
est_delay_by_AIF_correct_no_noisy_3D(CTC_4D_Mask_By_Val > 0) = est_delay_by_AIF_correct_3D(CTC_4D_Mask_By_Val > 0);

MeanFN=[Output_directory 'BAT_ACOPED_Noisy_' num2str(screen_val*100) '_Masked.nii'];
Raw2Nii(est_delay_by_AIF_correct_no_noisy_3D,MeanFN,'float32',DCEFNs{1});

est_delay_by_AIF_correct_no_noisy_brain_extract_3D      = zeros(size(est_delay_by_AIF_correct_no_noisy_3D));
est_delay_by_AIF_correct_no_noisy_brain_extract_3D(Brain_Mask_3D > 0) = est_delay_by_AIF_correct_no_noisy_3D(Brain_Mask_3D > 0);

MeanFN=[Output_directory 'BAT_ACOPED_Noisy_' num2str(screen_val*100) '_Masked_Brain_Extract.nii'];
Raw2Nii(est_delay_by_AIF_correct_no_noisy_brain_extract_3D,MeanFN,'float32',DCEFNs{1});


est_delay_by_AIF_correct_brain_extract_3D                    = zeros(size(est_delay_by_AIF_correct_3D));
est_delay_by_AIF_correct_brain_extract_3D(Brain_Mask_3D > 0) = est_delay_by_AIF_correct_3D(Brain_Mask_3D > 0);
MeanFN=[Output_directory 'BAT_ACOPED_Brain_Extract.nii'];
Raw2Nii(est_delay_by_AIF_correct_brain_extract_3D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Est_IRF_Larsson_Filter.nii'];
Raw2Nii(Est_IRF_no_Delay_4D,MeanFN,'float32',DCEFNs{1});

MeanFN=[Output_directory 'Conv_Res_with_IRF.nii'];
Raw2Nii(conv_result_no_Delay_IRF_4D,MeanFN,'float32',DCEFNs{1});

% Normalized 0-1 maps
Flow_no_Delay_3D          = loadniidata([Output_directory 'Flow_Larsson_no_Delay.nii']);
Ktrans_no_Delay_3D                = loadniidata([Output_directory 'Ktrans_no_Delay.nii']);
Vb_no_Delay_3D                    = loadniidata([Output_directory 'Vb_no_Delay.nii']);
Ve_no_Delay_3D                    = loadniidata([Output_directory 'Ve_no_Delay.nii']);

if Threshold_Norm_Maps
    max_val                                                            = Threshold_Val;
    Flow_Larsson_3D_Thresholded                                        = Flow_no_Delay_3D;
    Flow_Larsson_3D_Thresholded(Flow_Larsson_3D_Thresholded > max_val) = max_val;
    Flow_Larsson_3D_Norm_0_1                                           = Flow_Larsson_3D_Thresholded ./ max(max(max(Flow_Larsson_3D_Thresholded)));
else
    Flow_Larsson_3D_Norm_0_1                                           = Flow_no_Delay_3D ./ max(max(max(Flow_no_Delay_3D)));
end


Ktrans_3D_Norm_0_1       = Ktrans_no_Delay_3D ./ max(max(max(Ktrans_no_Delay_3D)));

MeanFN = [Output_directory 'Flow_Larsson_no_Delay_Normalized_0_1.nii'];
Raw2Nii(Flow_Larsson_3D_Norm_0_1,MeanFN,'float32',DCEFNs{1});

MeanFN = [Output_directory 'Ktrans_Larsson_no_Delay_Normalized_0_1.nii'];
Raw2Nii(Ktrans_3D_Norm_0_1,MeanFN,'float32',DCEFNs{1});


% Normalize maps (if masks exist)
if ( exist(WM_mask_absolute_path,'file') )
    
    display('-I- Normalizing Maps to White Matter...');
    
    WM_mask_3D                = loadniidata(WM_mask_absolute_path);
    
    % According to Larsson. WM Flow should be 30.6 [mL/100mL/min]
    [ Normalized_F_Map ]             = Normalize_Output_Maps( Flow_no_Delay_3D, WM_mask_3D , 30.6);
    
    % Take indices where both F and WM mask exist
    WM_mask_3D_Flow                           = zeros(size(WM_mask_3D));
    WM_mask_3D_Flow(F_Model_Selected_3D ~= 0) =  WM_mask_3D(F_Model_Selected_3D ~= 0);
    
    
    Normalized_F_Map_Brain_Extract                    = zeros(size(Normalized_F_Map));
    Normalized_F_Map_Brain_Extract(Brain_Mask_3D > 0) = Normalized_F_Map(Brain_Mask_3D > 0);
    
    MeanFN = [Output_directory 'Flow_Larsson_Relative_WM_30_6_Brain_Extract.nii'];
    Raw2Nii(Normalized_F_Map_Brain_Extract,MeanFN,'float32',DCEFNs{1});
    
    
    rCBF                = loadniidata('\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\forCorral\dce\rCBF.nii');
    rCBF_Brain_Extract                    = zeros(size(rCBF));
    rCBF_Brain_Extract(Brain_Mask_3D > 0) = rCBF(Brain_Mask_3D > 0);
    
    MeanFN = [Output_directory 'rCBF_Brain_Extract.nii'];
    Raw2Nii(rCBF_Brain_Extract,MeanFN,'float32',DCEFNs{1});
    
    
    
    [ Normalized_F_Model_Select_Map ]         = Normalize_Output_Maps( F_Model_Selected_3D, WM_mask_3D_Flow , 30.6);
    
    
    Normalized_F_Model_Select_Map_Brain_Extract  = zeros(size(Normalized_F_Model_Select_Map));
    Normalized_F_Model_Select_Map_Brain_Extract(Brain_Mask_3D > 0) = Normalized_F_Model_Select_Map(Brain_Mask_3D > 0);
    
    MeanFN = [Output_directory 'Flow_Larsson_model_selection_Relative_WM_30_6_Brain_Extract.nii'];
    Raw2Nii(Normalized_F_Model_Select_Map_Brain_Extract,MeanFN,'float32',DCEFNs{1});
    
    
    % According to Larsson. WM Ktrans should be 0.84 [mL/100mL/min]
    [ Normalized_Ktrans_Map ] = Normalize_Output_Maps( Ktrans_no_Delay_3D, WM_mask_3D , 0.84);
    % According to Jim.     WM Vp should be 0.01 [mL/100mL]
    [ Normalized_Vb_Map ]     = Normalize_Output_Maps( Vb_no_Delay_3D, WM_mask_3D , 0.01);
    
    MeanFN = [Output_directory 'Flow_Larsson_Relative_WM_30_6.nii'];
    Raw2Nii(Normalized_F_Map,MeanFN,'float32',DCEFNs{1});
    
    MeanFN = [Output_directory 'Flow_Larsson_model_selection_Relative_WM_30_6.nii'];
    Raw2Nii(Normalized_F_Model_Select_Map,MeanFN,'float32',DCEFNs{1});
    
    MeanFN = [Output_directory 'Ktrans_Relative_no_Delay_WM_0_84.nii'];
    Raw2Nii(Normalized_Ktrans_Map,MeanFN,'float32',DCEFNs{1});
    
    MeanFN = [Output_directory 'Vb_Relative_no_Delay_WM_0_01.nii'];
    Raw2Nii(Normalized_Vb_Map,MeanFN,'float32',DCEFNs{1});
    
end

end

