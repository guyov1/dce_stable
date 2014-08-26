% Close all open figures
close all;
% Clear everything except breakpoints
current_break_points = dbstatus;
save('myBreakpoints.mat', 'current_break_points');
clear all;
load('myBreakpoints.mat');
delete('myBreakpoints.mat');
dbstop(current_break_points);

% Add current directory and all it's subdirectories to Matlab's path
addpath(genpath('./'));

Verbosity = 'None';

% Initiate simulation struct
Sim_Struct = struct;

% Set simulation parameters
Sim_Struct = Simulation_Set_Params(Sim_Struct, Verbosity);

% Override the real data flag
Sim_Struct.RealData_Flag = true;

% Set parallel processing if needed
Set_Parallel_Processing(Sim_Struct, Verbosity);

% Read input MRI data
[Subject_name, Subject_Path, WM_mask_absolute_path, Art_Mask, After_CTC_mat, DCECoregP] = ReadRealData();

% Create output directory
mkdir([Subject_Path filesep 'Perfusion_DCE\']);
PefusionOutput        = [Subject_Path filesep 'Perfusion_DCE\'];
% Set output directory for figures/report
Output_directory      =  [PefusionOutput 'Run_Output/'];

% Add the needed data from the DCE run
display('-I- Loading .mat files from DCE run...');

if(exist(After_CTC_mat,'file'))
    load(After_CTC_mat);
else
    error('-E- AfterCTC.mat does not exist...');
end

% Define needed parameters from DCE data
Sim_Struct.num_time_stamps  = size(CTC2D,2);
Sim_Struct.num_voxels       = size(CTC2D,1);
Sim_Struct.sec_interval     = TimeBetweenDCEVolsFinal;
Sim_Struct.min_interval     = Sim_Struct.sec_interval / 60;
time_vec_minutes            = Sim_Struct.min_interval * ( 0 : Sim_Struct.num_time_stamps - 1 );
Sim_Struct.time_vec_minutes = time_vec_minutes;
 


%TimeBetweenDCEVolsMin   = TimeBetweenDCEVolsFinal/60;
%HSampleTs               = TimeBetweenDCEVolsMin*(0:nSVols-1);
%time_vec_seconds   = (1:nSVols).* sec_interval;
%time_vec_minutes   = time_vec_seconds / 60;
%nSVols                      = size(CTC2D,2);
%num_pixels                  = size(CTC2D,1);

%% Create AIF

if exist('AIFFindData_mat','var')
    load(AIFFindData_mat);
    
    % AIF_Parker8t=@(x,t) AIF_Parkerg2( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8))*x(2);
    AIF_Parker9t    =@(x,t)AIF_Parkerg3( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8),x(9),max(Sim_Struct.time_vec_minutes))*x(2);
    % Create Parker's AIF
    AIF_paramtertic =AIF_Parker9t(OutAIFParam,Sim_Struct.time_vec_minutes);
end


% Artreries and Veins masks

ICA_Art3D   = loadniidata(Art_Mask);
ICA_Art2D   = Reshape4d22d(ICA_Art3D,Msk2);
% Indices in 2D for Arteris/veins
ICA_Art2D_indices  = find(ICA_Art2D>0);

if exist('Vein_Mask','var')
    ICA_Vein3D  = loadniidata(Vein_Mask);
    ICA_Vein2D  = Reshape4d22d(ICA_Vein3D,Msk2);
    ICA_Vein2D_indices = find(ICA_Vein2D>0);
end

% Plot Arteris/veins Ct's and their averages
fig_num            = figure;

subplot(1,2,1);
if exist('AIFFindData_mat','var')
    plot(time_vec_minutes,AIF_paramtertic,time_vec_minutes,AIF_paramtertic,'*');
end
title('Functional AIF');
xlabel('Time [Min]');
ylabel('Amplitude');
subplot(1,2,2);
hold on;
plot(time_vec_minutes,transpose(CTC2D(ICA_Art2D_indices,:)),'b');
h1 = plot(time_vec_minutes,mean(transpose(CTC2D(ICA_Art2D_indices,:)),2),'k','LineWidth',4);
if exist('Vein_Mask','var')
    plot(time_vec_minutes,transpose(CTC2D(ICA_Vein2D_indices,:)),'g');
    h2 = plot(time_vec_minutes,mean(transpose(CTC2D(ICA_Vein2D_indices,:)),2),'r','LineWidth',4);
end

title('Arteris/Veins Cts and their averages');
xlabel('Time [Min]');
ylabel('Amplitude');
if exist('Vein_Mask','var')
    legend([h1 h2],'Mean Artery','Mean Vein');
else
    legend([h1],'Mean Artery');
end

hold off;

LogFN = [PefusionOutput 'Run_Output\Log.mat'];
if exist(LogFN,'file')
    delete(LogFN)
end
SN = 'Report';
Log.idx_000={['\\title{' SN '}\r\n\\maketitle\r\n']};
% Create directory if does not exist
if ~exist([PefusionOutput 'Run_Output'],'dir')
    mkdir([PefusionOutput 'Run_Output']);
end
save(LogFN,'Log');

gprint(fig_num,[PefusionOutput 'Run_Output/' 'InputAIFs.png']);
%gprint(fig_num,'Run_Output/InputAIFs.png');

AddToLog([PefusionOutput 'Run_Output\'],'idx_001','\\subsection*{Possible AIF Inputs}');
AddToLog([PefusionOutput 'Run_Output\'],'idx_002','InputAIFs','InputAIFs.png');

AIF_estimated_ICA = transpose(mean(transpose(CTC2D(ICA_Art2D_indices,:)),2));

%% Take Ct and AIF and calculate Ht

% intersect(find(sum(CTC2D,2)> 0.5),find(sum(CTC2D,2)< 0.6))
%Ct = CTC2D([44860 44865],:);
%Ct = CTC2D([40729],:); % Vein

%Ct = CTC2D([44865],:);
Ct = CTC2D(:,:);
DEBUG = true; % For drawing plots
%DEBUG = false; % For drawing plots

num_total_voxels = size(Ct,1);

% Choose the AIF (either parametric or from ICA average)
%Chosen_AIF = AIF_paramtertic;
Chosen_AIF = AIF_estimated_ICA;
%Chosen_AIF = transpose(smooth(AIF_estimated_ICA));

% Scale AIF as necessary
Chosen_AIF = Sim_Struct.AIF_Scaling_Factor * Chosen_AIF;

[ Flow_Larsson, Delay_sec_Larsson, t_delay_seconds, Delay_BiExp_Fit_seconds, sigma_seconds, amplitude, Est_ht, calculated_gaussian, conv_result_ht, conv_result_gaussian, RMS_ht, RMS_gauss, RMS_params,...
    calculated_double_gaussian, conv_result_double_gaussian, double_gauss_params, RMS_double_gauss, RMS_params_double_gauss, Ktrans, Vb, Ve, MTT, Ktrans_Patlak_vec, Vb_Patlak_vec, MTT_Patlak_vec ] = ...
    Get_Ht_Deconvolving(Sim_Struct, Chosen_AIF, Ct , Output_directory, Subject_name, Sim_Struct.Force_RealData_Calc, Verbosity);

%% Debugging results

if (DEBUG && num_total_voxels < 1000)
    
    fig_num = figure;
    
    % Plot estimated h(t) and calculated gaussian
    %hold on;plot(Est_ht,'b');plot(calculated_gaussian,'g');hold off;
    subplot(1,3,1);
    hold on;
    h1 = plot(time_vec_minutes,Chosen_AIF,'k');
    h2 = plot(time_vec_minutes,Chosen_AIF,'*');
    h3 = plot(time_vec_minutes,Ct,'b');
    h4 = plot(time_vec_minutes,Ct,'m*');
    hold off;
    title('Input AIF & Ct');
    xlabel('Time [Min]');
    legend([h2 h4],'Input AIF','Output Ct');
    
    subplot(1,3,2);
    hold on;
    h1 = plot(time_vec_minutes,Est_ht,'k');
    h2 = plot(time_vec_minutes,Est_ht,'g*');
    h3 = plot(time_vec_minutes,calculated_gaussian,'r');
    h4 = plot(time_vec_minutes,calculated_gaussian,'y*');
    h5 = plot(time_vec_minutes,calculated_double_gaussian,'m');
    h6 = plot(time_vec_minutes,calculated_double_gaussian,'c*');
    hold off;
    title('Estimated Filter h(t)');
    xlabel('Time [Min]')
    if ( size(h2,1)==1 )
        legend([h2 h4 h6],'Estimated h(t)','Fitted Gaussian', 'Fitted Double Gaussian');
    end
    
    subplot(1,3,3);
    hold on;
    h1 = plot(time_vec_minutes,Ct,'b');
    h2 = plot(time_vec_minutes,Ct,'m*');
    h3 = plot(time_vec_minutes,conv_result_ht,'k');
    h4 = plot(time_vec_minutes,conv_result_ht,'g*');
    h5 = plot(time_vec_minutes,conv_result_gaussian,'r');
    h6 = plot(time_vec_minutes,conv_result_gaussian,'y*');
    h7 = plot(time_vec_minutes,conv_result_double_gaussian,'m');
    h8 = plot(time_vec_minutes,conv_result_double_gaussian,'c*');
    hold off;
    title('Input Ct and est. Ct');
    xlabel('Time [Min]');
    if ( size(h2,1)==1 )
        legend([h2 h4 h6 h8],'Input Ct(t)','Estimated Ct(t) - h(t)','Estimated Ct(t) - Gaussian','Estimated Ct(t) - Double Gaussian');
    end
    
    gprint(fig_num,[PefusionOutput 'Run_Output/' 'VoxelAnalysis.png']);
    %gprint(fig_num,'Run_Output/VoxelAnalysis.png');
    
    AddToLog([PefusionOutput 'Run_Output\'],'idx_003','\\subsection*{Voxel Wise Analysis}');
    AddToLog([PefusionOutput 'Run_Output\'],'idx_004','VoxelAnalysis','VoxelAnalysis.png');
    
end

%% Save and Write results

display('-I- Saving parameters result of Get_Ht_Deconvolving...');

Mat_File_To_Save = [PefusionOutput 'Run_Output/' 'All_Parameters_Result.mat'];
%Mat_File_To_Save = 'Run_Output/All_Parameters_Result.mat';

save(Mat_File_To_Save,'Flow_Larsson','Delay_sec_Larsson'...
    ,'t_delay_seconds','Delay_BiExp_Fit_seconds','sigma_seconds','amplitude','Est_ht','conv_result_ht', 'conv_result_gaussian','RMS_ht'...
    ,'RMS_gauss','RMS_params','double_gauss_params','RMS_double_gauss','RMS_params_double_gauss','Ktrans', 'Vb', 'Ve', 'MTT', 'Ktrans_Patlak_vec', 'Vb_Patlak_vec', 'MTT_Patlak_vec' ,'Msk2'...
    ,'WorkingP','PefusionOutput','num_total_voxels','time_vec_minutes','TimeBetweenDCEVolsFinal','time_vec_minutes');
%load(Mat_File_To_Save);

% Write only in case we used the entire brain (say, above 1000 voxels)
if (num_total_voxels > 1000)
    
    % Reshape 2D to 3D and 4D
    t_delay_seconds_3D                   = Reshape2DCto4D(t_delay_seconds,Msk2);
    Delay_BiExp_Fit_seconds_3D           = Reshape2DCto4D(Delay_BiExp_Fit_seconds,Msk2);
    sigma_seconds_3D                     = Reshape2DCto4D(sigma_seconds,Msk2);
    amplitude_3D                         = Reshape2DCto4D(amplitude,Msk2);
    Est_ht_4D                            = Reshape2DCto4D(mat2cell(Est_ht,size(Est_ht,1),ones(1,size(Est_ht,2))),Msk2);
    CTC_4D                               = Reshape2DCto4D(mat2cell(CTC2D,size(CTC2D,1),ones(1,size(CTC2D,2))),Msk2);
    conv_result_ht_4D                    = Reshape2DCto4D(mat2cell(conv_result_ht,size(conv_result_ht,1),ones(1,size(conv_result_ht,2))),Msk2);
    conv_result_gaussian_3D              = Reshape2DCto4D(mat2cell(conv_result_gaussian,size(conv_result_gaussian,1),ones(1,size(conv_result_gaussian,2))),Msk2);
    RMS_ht_3D                            = Reshape2DCto4D(RMS_ht,Msk2);
    RMS_gauss_3D                         = Reshape2DCto4D(RMS_gauss,Msk2);
    RMS_params_3D                        = Reshape2DCto4D(RMS_params,Msk2);
    calculated_gaussian_3D               = Reshape2DCto4D(mat2cell(calculated_gaussian,size(calculated_gaussian,1),ones(1,size(calculated_gaussian,2))),Msk2);
    calculated_double_gaussian_3D        = Reshape2DCto4D(mat2cell(calculated_double_gaussian,size(calculated_double_gaussian,1),ones(1,size(calculated_double_gaussian,2))),Msk2);
    conv_result_double_gaussian_3D       = Reshape2DCto4D(mat2cell(conv_result_double_gaussian,size(conv_result_double_gaussian,1),ones(1,size(conv_result_double_gaussian,2))),Msk2);
    % Larsson's
    Flow_Larsson_3D                      = Reshape2DCto4D(Flow_Larsson,Msk2);
    Delay_sec_Larsson_3D                 = Reshape2DCto4D(Delay_sec_Larsson,Msk2);
    % Used AIF duplicated to all pixels
    AIF_Larsson_Duplicated_3D            = repmat(Chosen_AIF,max(size(Est_ht)),1);
    AIF_Larsson_4D                       = Reshape2DCto4D(mat2cell(AIF_Larsson_Duplicated_3D,size(AIF_Larsson_Duplicated_3D,1),ones(1,size(AIF_Larsson_Duplicated_3D,2))),Msk2);
    
    % Reshape double gaussian results
    t_delay_1_double_gauss_seconds_3D    = Reshape2DCto4D(double_gauss_params(1,:),Msk2);
    sigma_1_double_gauss_seconds_3D      = Reshape2DCto4D(double_gauss_params(2,:),Msk2);
    amplitude_1_double_gauss_seconds_3D  = Reshape2DCto4D(double_gauss_params(3,:),Msk2);
    t_delay_2_double_gauss_seconds_3D    = Reshape2DCto4D(double_gauss_params(4,:),Msk2);
    sigma_2_double_gauss_seconds_3D      = Reshape2DCto4D(double_gauss_params(5,:),Msk2);
    amplitude_2_double_gauss_seconds_3D  = Reshape2DCto4D(double_gauss_params(6,:),Msk2);
    RMS_double_gauss_3D                  = Reshape2DCto4D(RMS_double_gauss,Msk2);
    RMS_params_double_gauss_3D           = Reshape2DCto4D(RMS_params_double_gauss,Msk2);
    Ktrans_3D                                = Reshape2DCto4D(Ktrans,Msk2);
    Vb_3D                                = Reshape2DCto4D(Vb,Msk2);
    Ve_3D                                = Reshape2DCto4D(Ve,Msk2);
    MTT_3D                               = Reshape2DCto4D(MTT,Msk2);
    MTT_Patlak_3D                        = Reshape2DCto4D(MTT_Patlak_vec,Msk2);
    Ktrans_Patlak_3D                         = Reshape2DCto4D(Ktrans_Patlak_vec,Msk2);
    Vb_Patlak_3D                         = Reshape2DCto4D(Vb_Patlak_vec,Msk2);
    
    % Title for PDF before displaying images
    AddToLog([PefusionOutput 'Run_Output\'],'idx_005','\\subsection*{Parameters Maps}');
    
    % Reshape 2D to 4D
    fig_num = figure;
    subplot(1,2,1);
    imagesc(mritransform(t_delay_seconds_3D(:,:,2)));colorbar;%colormap(jet)
    title('Slice #2 - Time Delay [Sec]');
    subplot(1,2,2);
    imagesc(mritransform(t_delay_seconds_3D(:,:,3)));colorbar;
    title('Slice #3 - Time Delay [Sec]');
    gprint(fig_num,[PefusionOutput 'Run_Output/' 'Time_Delay_Gaussian.png']);
    %gprint(fig_num,'Run_Output/Time_Delay_Gaussian.png');
    AddToLog([PefusionOutput 'Run_Output\'],'idx_006','TimeDelayGaussian','Time_Delay_Gaussian.png');
    
    fig_num = figure;
    subplot(1,2,1);
    imagesc(mritransform(sigma_seconds_3D(:,:,2)));colorbar;
    title('Slice #2 - Sigma [Sec]');
    subplot(1,2,2);
    imagesc(mritransform(sigma_seconds_3D(:,:,3)));colorbar;
    title('Slice #3 - Sigma [Sec]');
    gprint(fig_num,[PefusionOutput 'Run_Output/' 'Sigma_Gaussian.png']);
    %gprint(fig_num,'Run_Output/Sigma_Gaussian.png');
    AddToLog([PefusionOutput 'Run_Output\'],'idx_007','SigmaGaussian','Sigma_Gaussian.png');
    
    fig_num = figure;
    subplot(1,2,1);
    imagesc(mritransform(RMS_params_3D(:,:,2)));colorbar;
    title('Slice #2 - RMS gaussian and ht');
    subplot(1,2,2);
    imagesc(mritransform(RMS_params_3D(:,:,3)));colorbar;
    title('Slice #3 - RMS gaussian and ht');
    gprint(fig_num,[PefusionOutput 'Run_Output/' 'RMS_gaussian_and_ht.png']);
    %gprint(fig_num,'Run_Output/RMS_gaussian_and_ht.png');
    AddToLog([PefusionOutput 'Run_Output\'],'idx_008','RmsGaussianAndHt','RMS_gaussian_and_ht.png');
    
    fig_num = figure;
    subplot(1,2,1);
    imagesc(mritransform(RMS_gauss_3D(:,:,2)));colorbar;
    title('Slice #2 - RMS gauss to data');
    subplot(1,2,2);
    imagesc(mritransform(RMS_gauss_3D(:,:,3)));colorbar;
    title('Slice #3 - RMS gauss to data');
    gprint(fig_num,[PefusionOutput 'Run_Output/' 'RMS_gauss_to_data.png']);
    %gprint(fig_num,'Run_Output/RMS_gauss_to_data.png');
    AddToLog([PefusionOutput 'Run_Output\'],'idx_009','RmsGaussToData','RMS_gauss_to_data.png');
    
    % Larsson's
    fig_num = figure;
    subplot(1,2,1);
    imagesc(mritransform(Flow_Larsson_3D(:,:,2)));colorbar;%colormap(jet)
    title('Slice #2 - Flow - Larson [mL/100g/min]');
    subplot(1,2,2);
    imagesc(mritransform(Flow_Larsson_3D(:,:,3)));colorbar;
    title('Slice #3 - Flow - Larson [mL/100g/min]');
    gprint(fig_num,[PefusionOutput 'Run_Output/' 'Flow_Larsson.png']);
    %gprint(fig_num,'Run_Output/Flow_Larsson.png');
    AddToLog([PefusionOutput 'Run_Output\'],'idx_010','FlowLarsson','Flow_Larsson.png');
    
    fig_num = figure;
    subplot(1,2,1);
    imagesc(mritransform(Delay_sec_Larsson_3D(:,:,2)));colorbar;%colormap(jet)
    title('Slice #2 - Delay - Larson [Sec]');
    subplot(1,2,2);
    imagesc(mritransform(Delay_sec_Larsson_3D(:,:,3)));colorbar;
    title('Slice #3 - Delay - Larson [Sec]');
    gprint(fig_num,[PefusionOutput 'Run_Output/' 'Delay_Larsson.png']);
    %gprint(fig_num,'Run_Output/Delay_Larsson.png');
    AddToLog([PefusionOutput 'Run_Output\'],'idx_011','DelayLarsson','Delay_Larsson.png');
    
    
    DDCE      = dir([DCECoregP '*.nii']);
    % HDR File
    DCEFNs    = strcat(DCECoregP,{DDCE.name})';
    
    MeanFN=[PefusionOutput 'FlowExtract_Time_Delay.nii'];
    Raw2Nii(t_delay_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Time_Delay_BiExp_Fit.nii'];
    Raw2Nii(Delay_BiExp_Fit_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Sigma.nii'];
    Raw2Nii(sigma_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Gaussian_Amplitude.nii'];
    Raw2Nii(amplitude_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Est_Ht_Larsson_Filter.nii'];
    Raw2Nii(Est_ht_4D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Conv_Res_with_Ht.nii'];
    Raw2Nii(conv_result_ht_4D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Est_Gaussian_Filter.nii'];
    Raw2Nii(calculated_gaussian_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Conv_Res_with_Gaussian.nii'];
    Raw2Nii(conv_result_gaussian_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Est_Double_Gaussian_Filter.nii'];
    Raw2Nii(calculated_double_gaussian_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Conv_Res_with_Double_Gaussian.nii'];
    Raw2Nii(conv_result_double_gaussian_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_RMS_Conv_Res_with_Ht_and_Ct.nii'];
    Raw2Nii(RMS_ht_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_RMS_Res_with_Gauss_and_Ct.nii'];
    Raw2Nii(RMS_gauss_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_RMS_Ht_with_Fitted_Gaussian.nii'];
    Raw2Nii(RMS_params_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_RMS_Res_with_Double_Gauss_and_Ct.nii'];
    Raw2Nii(RMS_double_gauss_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_RMS_Ht_with_Fitted_Double_Gaussian.nii'];
    Raw2Nii(RMS_params_double_gauss_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Ktrans.nii'];
    Raw2Nii(Ktrans_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Vb.nii'];
    Raw2Nii(Vb_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Ve.nii'];
    Raw2Nii(Ve_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_MTT.nii'];
    Raw2Nii(MTT_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_MTT_Patlak.nii'];
    Raw2Nii(MTT_Patlak_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Ktrans_Patlak.nii'];
    Raw2Nii(Ktrans_Patlak_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Vb_Patlak.nii'];
    Raw2Nii(Vb_Patlak_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Time_Delay_1_double_gaussian.nii'];
    Raw2Nii(t_delay_1_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Sigma_1_double_gaussian.nii'];
    Raw2Nii(sigma_1_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Gaussian_Amplitude_1_double_gaussian.nii'];
    Raw2Nii(amplitude_1_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Time_Delay_2_double_gaussian.nii'];
    Raw2Nii(t_delay_2_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Sigma_2_double_gaussian.nii'];
    Raw2Nii(sigma_2_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Gaussian_Amplitude_2_double_gaussian.nii'];
    Raw2Nii(amplitude_2_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Flow_Larsson.nii'];
    Raw2Nii(Flow_Larsson_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_Delay_Larsson.nii'];
    Raw2Nii(Delay_sec_Larsson_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_AIF_Used_Larsson.nii'];
    Raw2Nii(AIF_Larsson_4D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'FlowExtract_CTC_per_voxel.nii'];
    Raw2Nii(CTC_4D,MeanFN,'float32',DCEFNs{1});
end

display('-I- Normalizing Maps to White Matter...');
% Normalize maps (if masks exist)
if ( exist(WM_mask_absolute_path,'file') )
    
    WM_mask_3D = loadniidata(WM_mask_absolute_path);
    
    % DEBUG
    Flow_Larsson_3D = loadniidata([PefusionOutput 'FlowExtract_Flow_Larsson.nii']);
    Ktrans_3D           = loadniidata([PefusionOutput 'FlowExtract_Ktrans.nii']);
    
    if exist([WorkingP 'ManualArtNoBAT3\'],'dir')
        Ktrans_3D       = loadniidata([WorkingP 'ManualArtNoBAT3\KtransFinalN.nii']);
    else % Auto art
        Ktrans_3D       = loadniidata([WorkingP 'AutoArtBAT\KtransFinalN.nii']);
    end
    
    
    % According to Larsson. WM Flow should be 30.6 [mL/100mL]
    [ Normalized_F_Map ] = Normalize_Output_Maps( Flow_Larsson_3D, WM_mask_3D , 30.6);
    % According to Larsson. WM Ktrans should be 0.84 [mL/100mL/min]
    [ Normalized_Ktrans_Map ] = Normalize_Output_Maps( Ktrans_3D, WM_mask_3D , 0.84);
    
    MeanFN = [PefusionOutput 'FlowExtract_Flow_Larsson_Relative_WM_30_6.nii'];
    Raw2Nii(Normalized_F_Map,MeanFN,'float32',DCEFNs{1});
    MeanFN = [PefusionOutput 'FlowExtract_Ktrans_Relative_WM_0_84.nii'];
    Raw2Nii(Ktrans_3D,MeanFN,'float32',DCEFNs{1});
    
    % 0-1 maps
    max_val = 0.026768 * 242.9221;
    
    Flow_Larsson_3D_Thresholded = Flow_Larsson_3D;
    Flow_Larsson_3D_Thresholded(Flow_Larsson_3D_Thresholded > max_val) = max_val;
    Flow_Larsson_3D_Norm_0_1 = Flow_Larsson_3D_Thresholded ./ max(max(max(Flow_Larsson_3D_Thresholded)));
    
    Ktrans_3D_Norm_0_1           = Ktrans_3D ./ max(max(max(Ktrans_3D)));
    Ktrans_3D_Norm_0_1       = Ktrans_3D ./ max(max(max(Ktrans_3D)));
    
    MeanFN = [PefusionOutput 'FlowExtract_Flow_Larsson_Normalized_0_1.nii'];
    Raw2Nii(Flow_Larsson_3D_Norm_0_1,MeanFN,'float32',DCEFNs{1});
    MeanFN = [PefusionOutput 'FlowExtract_Ktrans_Larsson_Normalized_0_1.nii'];
    Raw2Nii(Ktrans_3D_Norm_0_1,MeanFN,'float32',DCEFNs{1});
    MeanFN = [PefusionOutput 'FlowExtract_Ktrans_Tofts_Normalized_0_1.nii'];
    Raw2Nii(Ktrans_3D_Norm_0_1,MeanFN,'float32',DCEFNs{1});
    
end

% Create PDF Report
Local_Path = [PefusionOutput 'Run_Output'];
MakeReport_func(Local_Path, LogFN);
%MakeReport;

% Display Ct(t) and fit
figure;
x = 187;
y = 149;
z = 2;
hold on;
plot(squeeze(CTC_4D(x,y,z,:)),'b');
plot(squeeze(conv_result_ht_4D(x,y,z,:)),'g');
hold off;


figure;
hold on;
h1 = plot(time_vec_minutes,squeeze(CTC_4D(x,y,z,:)),'LineWidth',6,'Color','k');
%h2 = plot(time_vec_minutes,AIF_part,'LineWidth',1,'LineStyle','+','Color','r');
%h3 = plot(time_vec_minutes,Kep_Filter_Part,'LineWidth',1,'LineStyle','o','Color','g');
%h4 = plot(time_vec_minutes,Sum_Result,'LineWidth',2,'Color','b');
%h5 = plot(time_vec_minutes,Sum_Result_NonLinear,'LineWidth',2,'Color','c');
hold off;
title('Concentration Time Curve','fontsize',15,'FontWeight','bold');
%legend([h1 h2 h3 h4],'CTC','Tofts AIF','Tofts permeability','Tofts fit')
xlabel('Time [Min]','fontsize',15,'FontWeight','bold');
ylabel('C_t(t) [mM]','fontsize',15,'FontWeight','bold');
set(gca,'fontsize',15,'FontWeight','bold');