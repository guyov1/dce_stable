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

% Override real data flag and parameters range
Sim_Struct.RealData_Flag       = true;
Sim_Struct.Vb_low              = 0.1; % [mL/100g]     , they used 3,6,12,18
Sim_Struct.Vb_max              = 100;
Sim_Struct.Ve_low              = 0.1; % Must be smaller than Vtis
Sim_Struct.Ve_max              = 100;
Sim_Struct.LowerBound_Larsson  = [Sim_Struct.Vb_low Sim_Struct.E_low  Sim_Struct.Ve_low];
Sim_Struct.UpperBound_Larsson  = [Sim_Struct.Vb_max Sim_Struct.E_max  Sim_Struct.Ve_max];
Sim_Struct.init_Ve_guess       = 0.1;
Sim_Struct.LQ_Model_AIF_Delay_Correct             = false;

% Set parallel processing if needed
Set_Parallel_Processing(Sim_Struct, Verbosity);

% Read input MRI data
[Subject_name, Subject_Path, WM_mask_absolute_path, Art_Mask, Vein_Mask, After_CTC_mat, DCECoregP] = ReadRealData();

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
    AIF_Parker9t    = @(x,t)AIF_Parkerg3( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8),x(9),max(Sim_Struct.time_vec_minutes))*x(2);
    % Create Parker's AIF
    AIF_paramtertic = AIF_Parker9t(OutAIFParam,Sim_Struct.time_vec_minutes);
end

% Artreries and Veins masks

ICA_Art3D   = loadniidata(Art_Mask);
ICA_Art2D   = Reshape4d22d(ICA_Art3D,Msk2);
% Indices in 2D for Arteris/veins
ICA_Art2D_indices  = find(ICA_Art2D>0);

if exist('Vein_Mask','var')
    ICA_Vein3D         = loadniidata(Vein_Mask);
    ICA_Vein2D         = Reshape4d22d(ICA_Vein3D,Msk2);
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

AIF_estimated_ICA  = transpose(mean(transpose(CTC2D(ICA_Art2D_indices,:)),2));

if exist('Vein_Mask','var')
    Vein_estimated_ICA = transpose(mean(transpose(CTC2D(ICA_Vein2D_indices,:)),2));
    % Coorect AIF scale to avoid Partial Volume Effect
    if Sim_Struct.Correct_PVE
        [AIF_estimated_ICA, Scale_Factor] = CorrectPVE(AIF_estimated_ICA, Vein_estimated_ICA, Sim_Struct);
    end
end

%% Take Ct and AIF and calculate Ht

% intersect(find(sum(CTC2D,2)> 0.5),find(sum(CTC2D,2)< 0.6))
%Ct = CTC2D([44860 44865],:);
%Ct = CTC2D([40729],:); % Vein

%Ct = CTC2D([44865],:);
Ct    = CTC2D(:,:);
DEBUG = true; % For drawing plots
%DEBUG = false; % For drawing plots

num_total_voxels = size(Ct,1);

% Choose the AIF (either parametric or from ICA average)
%Chosen_AIF = AIF_paramtertic;
Chosen_AIF = AIF_estimated_ICA;
%Chosen_AIF = transpose(smooth(AIF_estimated_ICA));

% Scale AIF as necessary
Chosen_AIF = double( Sim_Struct.AIF_Scaling_Factor * Chosen_AIF );

% TODO - add:
%               Vb_with_Delay_High_F
%               Vb_with_Delay_no_E
%               Vb_with_Delay_no_E_High_F
%
%               Ve_with_Delay_High_F
%               E_with_Delay_High_F
%
%               RMS_Larss_with_Delay_High_F
%               RMS_Larss_with_Delay_no_E_3D
%               RMS_Larss_with_Delay_no_E_High_F
%


[ Flow_Larsson_with_Delay, Flow_Larsson_no_Delay, Delay_sec_by_Max_Val_with_Delay, Delay_sec_by_Max_Val_no_Delay, est_delay_by_AIF_correct, t_delay_single_gauss_sec, ...
    sigma_seconds_single_gauss, Amp_single_gauss, Est_IRF_with_Delay, Est_IRF_no_Delay, fitted_gaussian, ...
    conv_result_Larss_with_Delay, conv_result_Larss_no_Delay, conv_result_Larss_no_Delay_High_F, conv_result_Larss_no_Delay_no_E,conv_result_Larss_no_Delay_no_E_High_F, ...
    conv_result_no_Delay_IRF, conv_result_gaussian, ...
    RMS_Larss_with_Delay, RMS_Larss_no_Delay, RMS_Larss_with_Delay_High_F, RMS_Larss_no_Delay_High_F, RMS_Larss_with_Delay_no_E, ...
    RMS_Larss_no_Delay_no_E, RMS_Larss_with_Delay_no_E_High_F, RMS_Larss_no_Delay_no_E_High_F, RMS_Larss_no_Delay_zero_params,...
    RMS_ht_no_Delay, RMS_gauss, RMS_params_Gauss, fitted_double_gaussian, conv_result_double_gaussian, ...
    double_gauss_params, RMS_double_gauss, RMS_params_double_gauss, Ktrans_with_Delay, Ktrans_no_Delay, ...
    E_with_Delay, E_no_Delay, Ktrans_with_Delay_High_F, Ktrans_no_Delay_High_F, ...
    Vb_with_Delay, Vb_no_Delay, Vb_with_Delay_High_F, Vb_no_Delay_High_F, Vb_with_Delay_no_E, Vb_no_Delay_no_E, Vb_with_Delay_no_E_High_F, Vb_no_Delay_no_E_High_F, ...
    Ve_with_Delay, Ve_no_Delay, Ve_with_Delay_High_F, Ve_no_Delay_High_F, ...
    MTT_with_Delay, MTT_no_Delay, Ktrans_Patlak_with_Delay_vec, Ktrans_Patlak_no_Delay_vec,...
    Vb_Patlak_with_Delay_vec, Vb_Patlak_no_Delay_vec, MTT_Patlak_with_Delay_vec, MTT_Patlak_no_Delay_vec ] = ...
    Get_Ht_Deconvolving(Sim_Struct, Chosen_AIF, Ct , Output_directory, Subject_name, Sim_Struct.Force_RealData_Calc, Verbosity);

%% Debugging results

if (DEBUG && num_total_voxels < 1000)
    
    fig_num = figure;
    
    % Plot estimated h(t) and calculated gaussian
    %hold on;plot(Est_IRF,'b');plot(calculated_gaussian,'g');hold off;
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
    h1 = plot(time_vec_minutes,Est_IRF_no_Delay,'k');
    h2 = plot(time_vec_minutes,Est_IRF_no_Delay,'g*');
    h3 = plot(time_vec_minutes,fitted_gaussian,'r');
    h4 = plot(time_vec_minutes,fitted_gaussian,'y*');
    h5 = plot(time_vec_minutes,fitted_double_gaussian,'m');
    h6 = plot(time_vec_minutes,fitted_double_gaussian,'c*');
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
    h3 = plot(time_vec_minutes,conv_result_no_Delay_IRF,'k');
    h4 = plot(time_vec_minutes,conv_result_no_Delay_IRF,'g*');
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

save(Mat_File_To_Save,'Flow_Larsson_with_Delay','Flow_Larsson_no_Delay','Delay_sec_by_Max_Val_with_Delay', 'Delay_sec_by_Max_Val_no_Delay'...
    ,'est_delay_by_AIF_correct','t_delay_single_gauss_sec','sigma_seconds_single_gauss','Amp_single_gauss'...
    ,'Est_IRF_with_Delay','Est_IRF_no_Delay', 'conv_result_Larss_with_Delay', 'conv_result_Larss_no_Delay', 'conv_result_Larss_no_Delay_High_F'...
    , 'conv_result_Larss_no_Delay_no_E','conv_result_Larss_no_Delay_no_E_High_F','conv_result_no_Delay_IRF', 'conv_result_gaussian'...
    , 'RMS_Larss_with_Delay', 'RMS_Larss_no_Delay', 'RMS_Larss_with_Delay_High_F', 'RMS_Larss_no_Delay_High_F', 'RMS_Larss_no_Delay_no_E','RMS_Larss_with_Delay_no_E', 'RMS_Larss_with_Delay_no_E_High_F', 'RMS_Larss_no_Delay_no_E_High_F','RMS_ht_no_Delay'...
    ,'RMS_gauss','RMS_params_Gauss','double_gauss_params','RMS_double_gauss','RMS_params_double_gauss'...
    ,'Ktrans_with_Delay','Ktrans_no_Delay', 'E_with_Delay', 'E_no_Delay', 'Ktrans_with_Delay_High_F', 'Ktrans_no_Delay_High_F','Vb_with_Delay', 'Vb_no_Delay','Vb_with_Delay_High_F'...
    , 'Vb_no_Delay_High_F', 'Vb_no_Delay_no_E', 'Vb_with_Delay_no_E', 'Vb_with_Delay_no_E_High_F', 'Vb_no_Delay_no_E_High_F'...
    , 'Ve_with_Delay', 'Ve_no_Delay','Ve_with_Delay_High_F', 'Ve_no_Delay_High_F', 'MTT_with_Delay', 'MTT_no_Delay'...
    , 'Ktrans_Patlak_with_Delay_vec', 'Ktrans_Patlak_no_Delay_vec'...
    , 'Vb_Patlak_with_Delay_vec', 'Vb_Patlak_no_Delay_vec', 'MTT_Patlak_with_Delay_vec', 'MTT_Patlak_no_Delay_vec'...
    ,'Msk2','WorkingP','PefusionOutput','num_total_voxels','time_vec_minutes','TimeBetweenDCEVolsFinal','time_vec_minutes');
%load(Mat_File_To_Save);

% Write only in case we used the entire brain (say, above 1000 voxels)
if (num_total_voxels > 1000)
    
    % Reshape 2D to 3D and 4D
    est_delay_by_AIF_correct_3D          = Reshape2DCto4D(est_delay_by_AIF_correct,Msk2);
    t_delay_single_gauss_sec_3D          = Reshape2DCto4D(t_delay_single_gauss_sec,Msk2);
    sigma_seconds_single_gauss_3D        = Reshape2DCto4D(sigma_seconds_single_gauss,Msk2);
    Amp_single_gauss_3D                  = Reshape2DCto4D(Amp_single_gauss,Msk2);
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
    Flow_Larsson_with_Delay_3D           = Reshape2DCto4D(Flow_Larsson_with_Delay,Msk2);
    Flow_Larsson_no_Delay_3D             = Reshape2DCto4D(Flow_Larsson_no_Delay,Msk2);
    Delay_sec_by_Max_Val_with_Delay_3D   = Reshape2DCto4D(Delay_sec_by_Max_Val_with_Delay,Msk2);
    Delay_sec_by_Max_Val_no_Delay_3D     = Reshape2DCto4D(Delay_sec_by_Max_Val_no_Delay,Msk2);
    % Used AIF duplicated to all pixels
    AIF_Larsson_Duplicated_3D            = repmat(Chosen_AIF,max(size(Est_IRF_no_Delay)),1);
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
    
    Ktrans_with_Delay_3D                 = Reshape2DCto4D(Ktrans_with_Delay,Msk2);
    Ktrans_no_Delay_3D                   = Reshape2DCto4D(Ktrans_no_Delay,Msk2);
    Ktrans_with_Delay_High_F_3D          = Reshape2DCto4D(Ktrans_with_Delay_High_F,Msk2);
    Ktrans_no_Delay_High_F_3D            = Reshape2DCto4D(Ktrans_no_Delay_High_F,Msk2);
    
    E_with_Delay_3D                      = Reshape2DCto4D(E_with_Delay,Msk2);
    E_no_Delay_3D                        = Reshape2DCto4D(E_no_Delay,Msk2);
    
    Vb_with_Delay_3D                     = Reshape2DCto4D(Vb_with_Delay,Msk2);
    Vb_no_Delay_3D                       = Reshape2DCto4D(Vb_no_Delay,Msk2);
    Vb_no_Delay_High_F_3D                = Reshape2DCto4D(Vb_no_Delay_High_F,Msk2);
    Vb_with_Delay_High_F_3D              = Reshape2DCto4D(Vb_with_Delay_High_F,Msk2);
    Vb_no_Delay_no_E_3D                  = Reshape2DCto4D(Vb_no_Delay_no_E,Msk2);
    Vb_with_Delay_no_E_3D                = Reshape2DCto4D(Vb_with_Delay_no_E,Msk2);
    Vb_no_Delay_no_E_High_F_3D           = Reshape2DCto4D(Vb_no_Delay_no_E_High_F,Msk2);
    Vb_with_Delay_no_E_High_F_3D         = Reshape2DCto4D(Vb_with_Delay_no_E_High_F,Msk2);
    
    Ve_with_Delay_3D                     = Reshape2DCto4D(Ve_with_Delay,Msk2);
    Ve_no_Delay_3D                       = Reshape2DCto4D(Ve_no_Delay,Msk2);
    Ve_with_Delay_High_F_3D              = Reshape2DCto4D(Ve_with_Delay_High_F,Msk2);
    Ve_no_Delay_High_F_3D                = Reshape2DCto4D(Ve_no_Delay_High_F,Msk2);
    
    MTT_with_Delay_3D                    = Reshape2DCto4D(MTT_no_Delay,Msk2);
    MTT_no_Delay_3D                      = Reshape2DCto4D(MTT_no_Delay,Msk2);
    
    MTT_with_Delay_Patlak_3D             = Reshape2DCto4D(MTT_Patlak_with_Delay_vec,Msk2);
    MTT_no_Delay_Patlak_3D               = Reshape2DCto4D(MTT_Patlak_no_Delay_vec,Msk2);
    Ktrans_Patlak_with_Delay_3D          = Reshape2DCto4D(Ktrans_Patlak_with_Delay_vec,Msk2);
    Ktrans_Patlak_no_Delay_3D            = Reshape2DCto4D(Ktrans_Patlak_no_Delay_vec,Msk2);
    Vb_Patlak_with_Delay_3D              = Reshape2DCto4D(Vb_Patlak_with_Delay_vec,Msk2);
    Vb_Patlak_no_Delay_3D                = Reshape2DCto4D(Vb_Patlak_no_Delay_vec,Msk2);
    
    %% PDF
    
    % Title for PDF before displaying images
    AddToLog([PefusionOutput 'Run_Output\'],'idx_005','\\subsection*{Parameters Maps}');
    
    % Reshape 2D to 4D
    fig_num = figure;
    subplot(1,2,1);
    imagesc(mritransform(est_delay_by_AIF_correct_3D(:,:,2)));colorbar;%colormap(jet)
    title('Slice #2 - Time Delay [Sec]');
    subplot(1,2,2);
    imagesc(mritransform(est_delay_by_AIF_correct_3D(:,:,3)));colorbar;
    title('Slice #3 - Time Delay [Sec]');
    gprint(fig_num,[PefusionOutput 'Run_Output/' 'Time_Delay_Novel_Method.png']);
    %gprint(fig_num,'Run_Output/Time_Delay_Gaussian.png');
    AddToLog([PefusionOutput 'Run_Output\'],'idx_006','TimeDelayNovel','Time_Delay_Novel_Method.png');
    
    fig_num = figure;
    subplot(1,2,1);
    imagesc(mritransform(sigma_seconds_single_gauss_3D(:,:,2)));colorbar;
    title('Slice #2 - Sigma [Sec]');
    subplot(1,2,2);
    imagesc(mritransform(sigma_seconds_single_gauss_3D(:,:,3)));colorbar;
    title('Slice #3 - Sigma [Sec]');
    gprint(fig_num,[PefusionOutput 'Run_Output/' 'Sigma_Gaussian.png']);
    %gprint(fig_num,'Run_Output/Sigma_Gaussian.png');
    AddToLog([PefusionOutput 'Run_Output\'],'idx_007','SigmaGaussian','Sigma_Gaussian.png');
    
    fig_num = figure;
    subplot(1,2,1);
    imagesc(mritransform(RMS_params_Gauss_3D(:,:,2)));colorbar;
    title('Slice #2 - RMS gaussian and ht');
    subplot(1,2,2);
    imagesc(mritransform(RMS_params_Gauss_3D(:,:,3)));colorbar;
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
    imagesc(mritransform(Flow_Larsson_no_Delay_3D(:,:,2)));colorbar;%colormap(jet)
    title('Slice #2 - Flow - Larson [mL/100g/min]');
    subplot(1,2,2);
    imagesc(mritransform(Flow_Larsson_no_Delay_3D(:,:,3)));colorbar;
    title('Slice #3 - Flow - Larson [mL/100g/min]');
    gprint(fig_num,[PefusionOutput 'Run_Output/' 'Flow_Larsson.png']);
    %gprint(fig_num,'Run_Output/Flow_Larsson.png');
    AddToLog([PefusionOutput 'Run_Output\'],'idx_010','FlowLarsson','Flow_Larsson.png');
    
    fig_num = figure;
    subplot(1,2,1);
    imagesc(mritransform(Delay_sec_by_Max_Val_no_Delay_3D(:,:,2)));colorbar;%colormap(jet)
    title('Slice #2 - Delay by Max Val- Larson [Sec]');
    subplot(1,2,2);
    imagesc(mritransform(Delay_sec_by_Max_Val_no_Delay_3D(:,:,3)));colorbar;
    title('Slice #3 - Delay by Max Val - Larson [Sec]');
    gprint(fig_num,[PefusionOutput 'Run_Output/' 'Delay_Larsson.png']);
    %gprint(fig_num,'Run_Output/Delay_Larsson.png');
    AddToLog([PefusionOutput 'Run_Output\'],'idx_011','DelayLarsson','Delay_Larsson.png');
    
    DDCE      = dir([DCECoregP filesep '*.nii']);
    % HDR File
    DCEFNs    = strcat([DCECoregP filesep],{DDCE.name})';
    
    %% Model Selection
    
    % First maps should be Gilad's
    %
    %
    
    nDataPoints      = length(time_vec_minutes);
    
    if Sim_Struct.Ignore_Delay_Model_Selection
        NParams          = [0 1 2 3 4];
        num_maps         = length(NParams);
        size_3D          = size(RMS_Larss_no_Delay_3D);
        RMSmaps          = zeros([size_3D num_maps]);
        
        RMSmaps(:,:,:,5) = RMS_Larss_with_Delay_3D;             % 5 Params
        %RMSmaps(:,:,:,8) = RMS_Larss_no_Delay_3D;               % 4 Params
        RMSmaps(:,:,:,4) = RMS_Larss_with_Delay_High_F_3D;      % 4 Params -----
        %RMSmaps(:,:,:,6) = RMS_Larss_no_Delay_High_F_3D;        % 3 Params
        RMSmaps(:,:,:,3) = RMS_Larss_with_Delay_no_E_3D;        % 3 Params -----
        %RMSmaps(:,:,:,4) = RMS_Larss_no_Delay_no_E_3D;          % 2 Params
        RMSmaps(:,:,:,2) = RMS_Larss_with_Delay_no_E_High_F_3D; % 2 Params -----
        %RMSmaps(:,:,:,2) = RMS_Larss_no_Delay_no_E_High_F_3D;   % 1 Params
        RMSmaps(:,:,:,1) = RMS_Larss_no_Delay_zero_params_3D;   % 0 Params
        
        [ ChosenByAIC_3D ]                      = Model_Selection( nDataPoints, RMSmaps, NParams, Sim_Struct.Data_Weight, Sim_Struct.AIC_Correction);
        %[ ChosenByAIC_3D ]                      = Model_Selection( nDataPoints, RMSmaps, NParams, 0.1, Sim_Struct.AIC_Correction);
        
        % F
        F_Model_Selected_3D                     = zeros(size(Flow_Larsson_with_Delay_3D));
        F_Model_Selected_3D (ChosenByAIC_3D==5) = Flow_Larsson_with_Delay_3D   (ChosenByAIC_3D==5);
        F_Model_Selected_3D (ChosenByAIC_3D==4) = zeros_map                    (ChosenByAIC_3D==4); % Put zero although its Inf
        F_Model_Selected_3D (ChosenByAIC_3D==3) = Flow_Larsson_with_Delay_3D   (ChosenByAIC_3D==3);        
        F_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2); % Put zero although its Inf
        F_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
        % Vb
        Vb_Model_Selected_3D                    = zeros(size(Vb_no_Delay_3D));
        Vb_Model_Selected_3D(ChosenByAIC_3D==5) = Vb_with_Delay_3D             (ChosenByAIC_3D==5);
        Vb_Model_Selected_3D(ChosenByAIC_3D==4) = Vb_with_Delay_High_F_3D      (ChosenByAIC_3D==4);
        Vb_Model_Selected_3D(ChosenByAIC_3D==3) = Vb_with_Delay_no_E_3D        (ChosenByAIC_3D==3);
        Vb_Model_Selected_3D(ChosenByAIC_3D==2) = Vb_with_Delay_no_E_High_F_3D (ChosenByAIC_3D==2);
        Vb_Model_Selected_3D(ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
        % Ve
        Ve_Model_Selected_3D                    = zeros(size(Ve_no_Delay_3D));
        Ve_Model_Selected_3D(ChosenByAIC_3D==5) = Ve_with_Delay_3D             (ChosenByAIC_3D==5);
        Ve_Model_Selected_3D(ChosenByAIC_3D==4) = Ve_with_Delay_High_F_3D      (ChosenByAIC_3D==4);
        Ve_Model_Selected_3D(ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
        Ve_Model_Selected_3D(ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
        Ve_Model_Selected_3D(ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
        % E
        E_Model_Selected_3D                     = zeros(size(E_no_Delay_3D));
        E_Model_Selected_3D (ChosenByAIC_3D==5) = E_with_Delay_3D              (ChosenByAIC_3D==5);
        E_Model_Selected_3D (ChosenByAIC_3D==4) = zeros_map                    (ChosenByAIC_3D==4); % Ktrans is not E
        E_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
        E_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
        E_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
        % Ktrans
        Ktrans_Model_Selected_3D                     = zeros(size(Ktrans_no_Delay_3D));
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==5) = Ktrans_with_Delay_3D         (ChosenByAIC_3D==5);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==4) = Ktrans_with_Delay_High_F_3D  (ChosenByAIC_3D==4);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
        
    else % Include Delay
        
        NParams          = [0 1 2 2 3 3 4 4 5];
        num_maps         = length(NParams);
        size_3D          = size(RMS_Larss_no_Delay_3D);
        RMSmaps          = zeros([size_3D num_maps]);
        
        RMSmaps(:,:,:,9) = RMS_Larss_with_Delay_3D;             % 5 Params
        RMSmaps(:,:,:,8) = RMS_Larss_no_Delay_3D;               % 4 Params
        RMSmaps(:,:,:,7) = RMS_Larss_with_Delay_High_F_3D;      % 4 Params -----
        RMSmaps(:,:,:,6) = RMS_Larss_no_Delay_High_F_3D;        % 3 Params
        RMSmaps(:,:,:,5) = RMS_Larss_with_Delay_no_E_3D;        % 3 Params -----
        RMSmaps(:,:,:,4) = RMS_Larss_no_Delay_no_E_3D;          % 2 Params
        RMSmaps(:,:,:,3) = RMS_Larss_with_Delay_no_E_High_F_3D; % 2 Params -----
        RMSmaps(:,:,:,2) = RMS_Larss_no_Delay_no_E_High_F_3D;   % 1 Params
        RMSmaps(:,:,:,1) = RMS_Larss_no_Delay_zero_params_3D;   % 0 Params
        
        [ ChosenByAIC_3D ]                      = Model_Selection( nDataPoints, RMSmaps, NParams, Sim_Struct.Data_Weight, Sim_Struct.AIC_Correction);
        %[ ChosenByAIC_3D ]                      = Model_Selection( nDataPoints, RMSmaps, NParams, 0.1, Sim_Struct.AIC_Correction);
        
        % Delay
        zeros_map                                      = zeros(size(Flow_Larsson_with_Delay_3D));
        Inf_map                                        = exp(100) * ones(size(Flow_Larsson_with_Delay_3D));
        AIF_Delay_Model_Selected_3D                    = zeros(size(est_delay_by_AIF_correct_3D));
        AIF_Delay_Model_Selected_3D(ChosenByAIC_3D==9) = est_delay_by_AIF_correct_3D (ChosenByAIC_3D==9);
        AIF_Delay_Model_Selected_3D(ChosenByAIC_3D==7) = est_delay_by_AIF_correct_3D (ChosenByAIC_3D==7);
        AIF_Delay_Model_Selected_3D(ChosenByAIC_3D==5) = est_delay_by_AIF_correct_3D (ChosenByAIC_3D==5);
        AIF_Delay_Model_Selected_3D(ChosenByAIC_3D==3) = est_delay_by_AIF_correct_3D (ChosenByAIC_3D==3);
        
        % F
        F_Model_Selected_3D                     = zeros(size(Flow_Larsson_with_Delay_3D));
        F_Model_Selected_3D (ChosenByAIC_3D==9) = Flow_Larsson_with_Delay_3D   (ChosenByAIC_3D==9);
        F_Model_Selected_3D (ChosenByAIC_3D==8) = Flow_Larsson_no_Delay_3D     (ChosenByAIC_3D==8);
        F_Model_Selected_3D (ChosenByAIC_3D==7) = zeros_map                    (ChosenByAIC_3D==7); % Put zero although its Inf
        F_Model_Selected_3D (ChosenByAIC_3D==6) = zeros_map                    (ChosenByAIC_3D==6); % Put zero although its Inf
        F_Model_Selected_3D (ChosenByAIC_3D==5) = Flow_Larsson_with_Delay_3D   (ChosenByAIC_3D==5);
        F_Model_Selected_3D (ChosenByAIC_3D==4) = Flow_Larsson_no_Delay_3D     (ChosenByAIC_3D==4);
        F_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3); % Put zero although its Inf
        F_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2); % Put zero although its Inf
        F_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
        
        % Vb
        Vb_Model_Selected_3D                    = zeros(size(Vb_no_Delay_3D));
        Vb_Model_Selected_3D(ChosenByAIC_3D==9) = Vb_with_Delay_3D             (ChosenByAIC_3D==9);
        Vb_Model_Selected_3D(ChosenByAIC_3D==8) = Vb_no_Delay_3D               (ChosenByAIC_3D==8);
        Vb_Model_Selected_3D(ChosenByAIC_3D==7) = Vb_with_Delay_High_F_3D      (ChosenByAIC_3D==7);
        Vb_Model_Selected_3D(ChosenByAIC_3D==6) = Vb_no_Delay_High_F_3D        (ChosenByAIC_3D==6);
        Vb_Model_Selected_3D(ChosenByAIC_3D==5) = Vb_with_Delay_no_E_3D        (ChosenByAIC_3D==5);
        Vb_Model_Selected_3D(ChosenByAIC_3D==4) = Vb_no_Delay_no_E_3D          (ChosenByAIC_3D==4);
        Vb_Model_Selected_3D(ChosenByAIC_3D==3) = Vb_with_Delay_no_E_High_F_3D (ChosenByAIC_3D==3);
        Vb_Model_Selected_3D(ChosenByAIC_3D==2) = Vb_no_Delay_no_E_High_F_3D   (ChosenByAIC_3D==2);
        Vb_Model_Selected_3D(ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
        % Ve
        Ve_Model_Selected_3D                    = zeros(size(Ve_no_Delay_3D));
        Ve_Model_Selected_3D(ChosenByAIC_3D==9) = Ve_with_Delay_3D             (ChosenByAIC_3D==9);
        Ve_Model_Selected_3D(ChosenByAIC_3D==8) = Ve_no_Delay_3D               (ChosenByAIC_3D==8);
        Ve_Model_Selected_3D(ChosenByAIC_3D==7) = Ve_with_Delay_High_F_3D      (ChosenByAIC_3D==7);
        Ve_Model_Selected_3D(ChosenByAIC_3D==6) = Ve_no_Delay_High_F_3D        (ChosenByAIC_3D==6);
        Ve_Model_Selected_3D(ChosenByAIC_3D==5) = zeros_map                    (ChosenByAIC_3D==5);
        Ve_Model_Selected_3D(ChosenByAIC_3D==4) = zeros_map                    (ChosenByAIC_3D==4);
        Ve_Model_Selected_3D(ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
        Ve_Model_Selected_3D(ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
        Ve_Model_Selected_3D(ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
        % E
        E_Model_Selected_3D                     = zeros(size(E_no_Delay_3D));
        E_Model_Selected_3D (ChosenByAIC_3D==9) = E_with_Delay_3D              (ChosenByAIC_3D==9);
        E_Model_Selected_3D (ChosenByAIC_3D==8) = E_no_Delay_3D                (ChosenByAIC_3D==8);
        E_Model_Selected_3D (ChosenByAIC_3D==7) = zeros_map                    (ChosenByAIC_3D==7); % Ktrans is not E
        E_Model_Selected_3D (ChosenByAIC_3D==6) = zeros_map                    (ChosenByAIC_3D==6); % Ktrans is not E
        E_Model_Selected_3D (ChosenByAIC_3D==5) = zeros_map                    (ChosenByAIC_3D==5);
        E_Model_Selected_3D (ChosenByAIC_3D==4) = zeros_map                    (ChosenByAIC_3D==4);
        E_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
        E_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
        E_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
        
        % Ktrans
        Ktrans_Model_Selected_3D                     = zeros(size(Ktrans_no_Delay_3D));
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==9) = Ktrans_with_Delay_3D         (ChosenByAIC_3D==9);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==8) = Ktrans_no_Delay_3D           (ChosenByAIC_3D==8);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==7) = Ktrans_with_Delay_High_F_3D  (ChosenByAIC_3D==7);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==6) = Ktrans_no_Delay_High_F_3D    (ChosenByAIC_3D==6);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==5) = zeros_map                    (ChosenByAIC_3D==5);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==4) = zeros_map                    (ChosenByAIC_3D==4);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==3) = zeros_map                    (ChosenByAIC_3D==3);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==2) = zeros_map                    (ChosenByAIC_3D==2);
        Ktrans_Model_Selected_3D (ChosenByAIC_3D==1) = zeros_map                    (ChosenByAIC_3D==1);
        
    end
    
    % Model Selection Map
    MeanFN=[PefusionOutput 'Model_Selection_Map.nii'];
    Raw2Nii(ChosenByAIC_3D,MeanFN,'float32',DCEFNs{1});
    
    %% Maps
    
    if (Sim_Struct.USE_ONE_GAUSSIAN)
        MeanFN=[PefusionOutput 'Delay_Single_Gaussian.nii'];
        Raw2Nii(t_delay_single_gauss_sec_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'Sigma_Single_Gaussian.nii'];
        Raw2Nii(sigma_seconds_single_gauss_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'Single_Gaussian_Amplitude.nii'];
        Raw2Nii(Amp_single_gauss_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'Est_Single_Gaussian_Filter.nii'];
        Raw2Nii(calculated_gaussian_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'Conv_Res_with_Single_Gaussian.nii'];
        Raw2Nii(conv_result_gaussian_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'RMS_Res_with_Single_Gauss_and_Ct.nii'];
        Raw2Nii(RMS_gauss_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'RMS_Ht_with_Fitted_Single_Gaussian.nii'];
        Raw2Nii(RMS_params_Gauss_3D,MeanFN,'float32',DCEFNs{1});
        
    end
    
    if (Sim_Struct.USE_DOUBLE_GAUSSIAN)
        
        MeanFN=[PefusionOutput 'Time_Delay_1_double_gaussian.nii'];
        Raw2Nii(t_delay_1_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'Sigma_1_double_gaussian.nii'];
        Raw2Nii(sigma_1_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'Amplitude_1_double_gaussian.nii'];
        Raw2Nii(amplitude_1_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'Time_Delay_2_double_gaussian.nii'];
        Raw2Nii(t_delay_2_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'Sigma_2_double_gaussian.nii'];
        Raw2Nii(sigma_2_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'Amplitude_2_double_gaussian.nii'];
        Raw2Nii(amplitude_2_double_gauss_seconds_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'Est_Double_Gaussian_Filter.nii'];
        Raw2Nii(calculated_double_gaussian_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'Conv_Res_with_Double_Gaussian.nii'];
        Raw2Nii(conv_result_double_gaussian_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'RMS_Res_with_Double_Gauss_and_Ct.nii'];
        Raw2Nii(RMS_double_gauss_3D,MeanFN,'float32',DCEFNs{1});
        
        MeanFN=[PefusionOutput 'RMS_Ht_with_Fitted_Double_Gaussian.nii'];
        Raw2Nii(RMS_params_double_gauss_3D,MeanFN,'float32',DCEFNs{1});
        
    end
    
    % ------------------- Delay -----------------------------
    MeanFN=[PefusionOutput 'Time_Delay_Novel_AIF_Correct.nii'];
    Raw2Nii(est_delay_by_AIF_correct_3D,MeanFN,'float32',DCEFNs{1});
    
    if ~Sim_Struct.Ignore_Delay_Model_Selection
        MeanFN=[PefusionOutput 'Time_Delay_Novel_AIF_Correct_Model_Selection.nii'];
        Raw2Nii(AIF_Delay_Model_Selected_3D,MeanFN,'float32',DCEFNs{1});
    end
    
    MeanFN=[PefusionOutput 'Delay_by_Max_Val_Larsson_no_Delay.nii'];
    Raw2Nii(Delay_sec_by_Max_Val_no_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Delay_by_Max_Val_Larsson_with_Delay.nii'];
    Raw2Nii(Delay_sec_by_Max_Val_with_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    % -------------------- F ----------------------------
    
    MeanFN=[PefusionOutput 'Flow_Larsson_with_Delay.nii'];
    Raw2Nii(Flow_Larsson_with_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Flow_Larsson_no_Delay.nii'];
    Raw2Nii(Flow_Larsson_no_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Flow_Model_Selection.nii'];
    Raw2Nii(F_Model_Selected_3D,MeanFN,'float32',DCEFNs{1});
    
    % -------------------- E ----------------------------
    
    MeanFN=[PefusionOutput 'E_with_Delay.nii'];
    Raw2Nii(E_with_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'E_no_Delay.nii'];
    Raw2Nii(E_no_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    %MeanFN=[PefusionOutput 'E_with_Delay_High_F.nii'];
    %Raw2Nii(Ktrans_with_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});
    
    %MeanFN=[PefusionOutput 'E_no_Delay_High_F.nii'];
    %Raw2Nii(Ktrans_no_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'E_Model_Selected.nii'];
    Raw2Nii(E_Model_Selected_3D,MeanFN,'float32',DCEFNs{1});
    
    % -------------------- Ktrans ----------------------------
    
    MeanFN=[PefusionOutput 'Ktrans_with_Delay.nii'];
    Raw2Nii(Ktrans_with_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Ktrans_no_Delay.nii'];
    Raw2Nii(Ktrans_no_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Ktrans_with_Delay_High_F.nii'];
    Raw2Nii(Ktrans_with_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Ktrans_no_Delay_High_F.nii'];
    Raw2Nii(Ktrans_no_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Ktrans_with_Delay_Patlak.nii'];
    Raw2Nii(Ktrans_Patlak_with_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Ktrans_no_Delay_Patlak.nii'];
    Raw2Nii(Ktrans_Patlak_no_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Ktrans_Model_Selected.nii'];
    Raw2Nii(Ktrans_Model_Selected_3D,MeanFN,'float32',DCEFNs{1});
    
    % -------------------- Vb ----------------------------
    
    MeanFN=[PefusionOutput 'Vb_with_Delay.nii'];
    Raw2Nii(Vb_with_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Vb_no_Delay.nii'];
    Raw2Nii(Vb_no_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Vb_with_Delay_High_F.nii'];
    Raw2Nii(Vb_with_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Vb_no_Delay_High_F.nii'];
    Raw2Nii(Vb_no_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Vb_with_Delay_no_E.nii'];
    Raw2Nii(Vb_with_Delay_no_E_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Vb_no_Delay_no_E.nii'];
    Raw2Nii(Vb_no_Delay_no_E_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Vb_with_Delay_no_E_High_F.nii'];
    Raw2Nii(Vb_with_Delay_no_E_High_F_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Vb_no_Delay_no_E_High_F.nii'];
    Raw2Nii(Vb_no_Delay_no_E_High_F_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Vb_Model_Selection.nii'];
    Raw2Nii(Vb_Model_Selected_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Vb_with_Delay_Patlak.nii'];
    Raw2Nii(Vb_Patlak_with_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Vb_no_Delay_Patlak.nii'];
    Raw2Nii(Vb_Patlak_no_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    % -------------------- Ve ----------------------------
    
    MeanFN=[PefusionOutput 'Ve_with_Delay.nii'];
    Raw2Nii(Ve_with_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Ve_no_Delay.nii'];
    Raw2Nii(Ve_no_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Ve_with_Delay_High_F.nii'];
    Raw2Nii(Ve_with_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Ve_no_Delay_High_F.nii'];
    Raw2Nii(Ve_no_Delay_High_F_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Ve_Model_Selection.nii'];
    Raw2Nii(Ve_Model_Selected_3D,MeanFN,'float32',DCEFNs{1});
    
    % -------------------- MTT ----------------------------
    
    MeanFN=[PefusionOutput 'MTT_with_Delay.nii'];
    Raw2Nii(MTT_with_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'MTT_no_Delay.nii'];
    Raw2Nii(MTT_no_Delay_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'MTT_with_Delay_Patlak.nii'];
    Raw2Nii(MTT_with_Delay_Patlak_3D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'MTT_no_Delay_Patlak.nii'];
    Raw2Nii(MTT_no_Delay_Patlak_3D,MeanFN,'float32',DCEFNs{1});
    
    % -------------------- RMS ----------------------------
    
    MeanFN=[PefusionOutput 'Log_RMS_Conv_Res_with_IRF_and_Ct.nii'];
    Raw2Nii(log(RMS_ht_no_Delay_3D),MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Log_RMS_Larss_with_Delay.nii'];
    Raw2Nii(log(RMS_Larss_with_Delay_3D),MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Log_RMS_Larss_No_Delay.nii'];
    Raw2Nii(log(RMS_Larss_no_Delay_3D),MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Log_RMS_Larss_with_Delay_High_F.nii'];
    Raw2Nii(log(RMS_Larss_with_Delay_High_F_3D),MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Log_RMS_Larss_no_Delay_High_F.nii'];
    Raw2Nii(log(RMS_Larss_no_Delay_High_F_3D),MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Log_RMS_Larss_with_Delay_no_E.nii'];
    Raw2Nii(log(RMS_Larss_with_Delay_no_E_3D),MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Log_RMS_Larss_no_Delay_no_E.nii'];
    Raw2Nii(log(RMS_Larss_no_Delay_no_E_3D),MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Log_RMS_Larss_with_Delay_no_E_High_F.nii'];
    Raw2Nii(log(RMS_Larss_with_Delay_no_E_High_F_3D),MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Log_RMS_Larss_no_Delay_no_E_High_F.nii'];
    Raw2Nii(log(RMS_Larss_no_Delay_no_E_High_F_3D),MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Log_RMS_Larss_no_Delay_zero_params.nii'];
    Raw2Nii(log(RMS_Larss_no_Delay_zero_params_3D),MeanFN,'float32',DCEFNs{1});
    
    % -------------------- General ----------------------------
    
    MeanFN=[PefusionOutput 'AIF_Used_Larsson.nii'];
    Raw2Nii(AIF_Larsson_4D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'CTC_per_voxel.nii'];
    Raw2Nii(CTC_4D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Est_IRF_Larsson_Filter.nii'];
    Raw2Nii(Est_IRF_no_Delay_4D,MeanFN,'float32',DCEFNs{1});
    
    MeanFN=[PefusionOutput 'Conv_Res_with_IRF.nii'];
    Raw2Nii(conv_result_no_Delay_IRF_4D,MeanFN,'float32',DCEFNs{1});
    
    
    
    
    % Normalized 0-1 maps
    Flow_Larsson_no_Delay_3D          = loadniidata([PefusionOutput 'Flow_Larsson_no_Delay.nii']);
    Ktrans_no_Delay_3D                = loadniidata([PefusionOutput 'Ktrans_no_Delay.nii']);
    Vb_no_Delay_3D                    = loadniidata([PefusionOutput 'Vb_no_Delay.nii']);
    Ve_no_Delay_3D                    = loadniidata([PefusionOutput 'Ve_no_Delay.nii']);
    
    if Sim_Struct.Threshold_Norm_Maps
        max_val                                                            = Sim_Struct.Threshold_Val;
        Flow_Larsson_3D_Thresholded                                        = Flow_Larsson_no_Delay_3D;
        Flow_Larsson_3D_Thresholded(Flow_Larsson_3D_Thresholded > max_val) = max_val;
        Flow_Larsson_3D_Norm_0_1                                           = Flow_Larsson_3D_Thresholded ./ max(max(max(Flow_Larsson_3D_Thresholded)));
    else
        Flow_Larsson_3D_Norm_0_1                                           = Flow_Larsson_no_Delay_3D ./ max(max(max(Flow_Larsson_no_Delay_3D)));
    end
    
    
    Ktrans_3D_Norm_0_1       = Ktrans_no_Delay_3D ./ max(max(max(Ktrans_no_Delay_3D)));
    
    MeanFN = [PefusionOutput 'Flow_Larsson_no_Delay_Normalized_0_1.nii'];
    Raw2Nii(Flow_Larsson_3D_Norm_0_1,MeanFN,'float32',DCEFNs{1});
    
    MeanFN = [PefusionOutput 'Ktrans_Larsson_no_Delay_Normalized_0_1.nii'];
    Raw2Nii(Ktrans_3D_Norm_0_1,MeanFN,'float32',DCEFNs{1});
    
    
    % Normalize maps (if masks exist)
    if ( exist(WM_mask_absolute_path,'file') )
        
        display('-I- Normalizing Maps to White Matter...');
        
        WM_mask_3D                = loadniidata(WM_mask_absolute_path);
        
        % According to Larsson. WM Flow should be 30.6 [mL/100mL/min]
        [ Normalized_F_Map ]      = Normalize_Output_Maps( Flow_Larsson_no_Delay_3D, WM_mask_3D , 30.6);
        % According to Larsson. WM Ktrans should be 0.84 [mL/100mL/min]
        [ Normalized_Ktrans_Map ] = Normalize_Output_Maps( Ktrans_no_Delay_3D, WM_mask_3D , 0.84);
        % According to Jim.     WM Vp should be 0.01 [mL/100mL]
        [ Normalized_Vb_Map ]     = Normalize_Output_Maps( Vb_no_Delay_3D, WM_mask_3D , 0.01);
        
        MeanFN = [PefusionOutput 'Flow_Larsson_no_Delay_Relative_WM_30_6.nii'];
        Raw2Nii(Normalized_F_Map,MeanFN,'float32',DCEFNs{1});
        
        MeanFN = [PefusionOutput 'Ktrans_Relative_no_Delay_WM_0_84.nii'];
        Raw2Nii(Normalized_Ktrans_Map,MeanFN,'float32',DCEFNs{1});
        
        MeanFN = [PefusionOutput 'Vb_Relative_no_Delay_WM_0_01.nii'];
        Raw2Nii(Normalized_Vb_Map,MeanFN,'float32',DCEFNs{1});
        
    end
    
end

% Create PDF Report
Local_Path = [PefusionOutput 'Run_Output'];
MakeReport_func(Local_Path, LogFN);
close all;

%MakeReport;

% % Display Ct(t) and fit of a single voxel
% figure;
% x = 187;
% y = 149;
% z = 2;
% hold on;
% plot(squeeze(CTC_4D(x,y,z,:)),'b');
% plot(squeeze(conv_result_IRF_4D(x,y,z,:)),'g');
% hold off;
%
% figure;
% hold on;
% h1 = plot(time_vec_minutes,squeeze(CTC_4D(x,y,z,:)),'LineWidth',6,'Color','k');
% %h2 = plot(time_vec_minutes,AIF_part,'LineWidth',1,'LineStyle','+','Color','r');
% %h3 = plot(time_vec_minutes,Kep_Filter_Part,'LineWidth',1,'LineStyle','o','Color','g');
% %h4 = plot(time_vec_minutes,Sum_Result,'LineWidth',2,'Color','b');
% %h5 = plot(time_vec_minutes,Sum_Result_NonLinear,'LineWidth',2,'Color','c');
% hold off;
% title('Concentration Time Curve','fontsize',15,'FontWeight','bold');
% %legend([h1 h2 h3 h4],'CTC','Tofts AIF','Tofts permeability','Tofts fit')
% xlabel('Time [Min]','fontsize',15,'FontWeight','bold');
% ylabel('C_t(t) [mM]','fontsize',15,'FontWeight','bold');
% set(gca,'fontsize',15,'FontWeight','bold');
%
%
%
