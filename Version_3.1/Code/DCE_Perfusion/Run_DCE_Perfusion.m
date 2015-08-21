function [ ] = Run_DCE_Perfusion( Subject_name, Subject_Path, Sim_Struct, PefusionOutput, Verbosity )
%Run_DCE_Perfusion Run Full DCE Perfusion Analysis

display('---------------------------------------------------');
display('-I- Started Run_DCE_Perfusion execution at:');
c = clock;
disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
display('---------------------------------------------------');

DCECoregP             = Subject_Path;
WM_mask_absolute_path = [Subject_Path  '\RefT1_WM_830.nii'];
Art_Mask              = [Subject_Path  '\InspectedRepVox.nii'];
Vein_Mask             = [Subject_Path  '\Veins_Mask.nii'];
After_CTC_mat         = [Subject_Path  '\AfterCTC.mat'];
Brain_Extract_path    = [Subject_Path  '\Manual_BrainMask.nii'];
Manual_AIF            = [Subject_Path  '\AIF.CSV'];

% Brain Mask
Brain_Mask_3D         = loadniidata(Brain_Extract_path);

% Override real data flag and parameters range
Sim_Struct.RealData_Flag              = true;
Sim_Struct.Vb_low                     = 0.1; % [mL/100g]     , they used 3,6,12,18
Sim_Struct.Vb_max                     = 100;
Sim_Struct.Ve_low                     = 0.1; % Must be smaller than Vtis
Sim_Struct.Ve_max                     = 100;
Sim_Struct.LowerBound_Larsson         = [Sim_Struct.Vb_low Sim_Struct.E_low  Sim_Struct.Ve_low];
Sim_Struct.UpperBound_Larsson         = [Sim_Struct.Vb_max Sim_Struct.E_max  Sim_Struct.Ve_max];
Sim_Struct.init_Ve_guess              = 0.1;

% Set parallel processing if needed
if Sim_Struct.Parallel_Real_Data_Est
    Set_Parallel_Processing(Sim_Struct, Verbosity);
end

% Set output directory for figures/report
Output_directory      =  [PefusionOutput 'Run_Output' filesep];
% Create directory if does not exist
if ~exist(Output_directory,'dir')
    mkdir(Output_directory);
end

% Add the needed data from the DCE run
display('-I- Loading .mat files from DCE run...');

if(exist(After_CTC_mat,'file'))
    load(After_CTC_mat);
else
    error('-E- AfterCTC.mat does not exist...');
end

%% Create AIF
if Sim_Struct.manual_aif
    
    % When time vector is not in equal distances, interpolate
    
    % Read the AIF.csv file
    tmp_aif_matrix                = csvread(Manual_AIF);
    
    % Assign to AIF struct
    AIF_Struct                    = struct();
    AIF_Struct.AIF_estimated_ICA  = tmp_aif_matrix(1,:);
    AIF_Struct.AIF_parametric     = tmp_aif_matrix(1,:);
    
    % Get time vector
    time_vec_minutes              = tmp_aif_matrix(2,:);
    time_diff_vec                 = diff(time_vec_minutes);
    
    % Check the diff vector is bigger than zero
    if min(time_diff_vec) < 0 
        error('-E- Error in input time vector of manual AIF!');
    end
    
    % Interpolate only in case the time difference vector is not equal
    if ~all(time_diff_vec == time_diff_vec(1))
        
        % New interval will be the biggest possible for equal distances
        new_time_interval = min(time_diff_vec);
        first_time_point  = time_vec_minutes(1);
        last_time_point   = time_vec_minutes(end);
        
        % Create a new time vector
        new_time_vector   = first_time_point : new_time_interval : last_time_point;
        
        % Interpolate old data using shape-preserving piecewise cubic
        % interpolation (each CTC should be in each vector)
        display('-I- Interpolating CTC...');
        new_CTC                     = interp1(time_vec_minutes, CTC2D', new_time_vector, 'PCHIP');
        
        % Update the old vectors and matrices
        Sim_Struct.min_interval     = new_time_interval;
        Sim_Struct.sec_interval     = new_time_interval * 60;
        Sim_Struct.time_vec_minutes = new_time_vector;
        
        % Overwrite CTC2D
        CTC2D                       = new_CTC';
        
    else
        Sim_Strct.sec_interval      = TimeBetweenDCEVolsFinal;
        Sim_Struct.min_interval     = Sim_Struct.sec_interval / 60;
        time_vec_minutes            = Sim_Struct.min_interval * ( 0 : Sim_Struct.num_time_stamps - 1 );
        Sim_Struct.time_vec_minutes = time_vec_minutes;
    end

    % Use interpolated data
    Sim_Struct.num_time_stamps  = size(CTC2D,2);
    Sim_Struct.num_voxels       = size(CTC2D,1);
else
    % Define needed parameters from DCE data
    Sim_Struct.num_time_stamps  = size(CTC2D,2);
    Sim_Struct.num_voxels       = size(CTC2D,1);
    Sim_Struct.sec_interval     = TimeBetweenDCEVolsFinal;
    Sim_Struct.min_interval     = Sim_Struct.sec_interval / 60;
    time_vec_minutes            = Sim_Struct.min_interval * ( 0 : Sim_Struct.num_time_stamps - 1 );
    Sim_Struct.time_vec_minutes = time_vec_minutes;

    if exist('AIFFindData_mat','var')
        [AIF_Struct] = chooseAifForRealData(Sim_Struct, CTC2D, Art_Mask, Vein_Mask, Msk2, Output_directory, AIFFindData_mat);
    else
        [AIF_Struct] = chooseAifForRealData(Sim_Struct, CTC2D, Art_Mask, Vein_Mask, Msk2, Output_directory);
    end
end

%% Take Ct and AIF and calculate Ht
Ct               = CTC2D(:,:);
num_total_voxels = Sim_Struct.num_voxels;

% Choose the AIF (either parametric or from ICA average)
Chosen_AIF = AIF_Struct.AIF_estimated_ICA; % AIF_paramtertic, transpose(smooth(AIF_estimated_ICA))

% Scale AIF as necessary
Chosen_AIF = double( Sim_Struct.AIF_Scaling_Factor * Chosen_AIF );

[ resultStruct ] = Get_Ht_Deconvolving(Sim_Struct, Chosen_AIF, Ct , Output_directory, Subject_name, Sim_Struct.Force_RealData_Calc, Verbosity);

%% Save and Write results

display('-I- Saving parameters result of Get_Ht_Deconvolving...');

Mat_File_To_Save = [Output_directory 'All_Parameters_Result.mat'];

save(Mat_File_To_Save,'resultStruct','Msk2','WorkingP','PefusionOutput','num_total_voxels','time_vec_minutes','TimeBetweenDCEVolsFinal','time_vec_minutes'...
    ,'Chosen_AIF', 'DCECoregP','Sim_Struct','WM_mask_absolute_path','Subject_Path');

% Write only in case we used the entire brain (say, above 1000 voxels)
if (num_total_voxels > 1000)
    resultToNiiAndReport(resultStruct, time_vec_minutes, CTC2D, Chosen_AIF, Msk2, Brain_Mask_3D, Output_directory, DCECoregP, Sim_Struct, WM_mask_absolute_path);
end

%% Create PDF Report
LogFN = [WorkingP 'Log.mat'];
MakeReport_func(Output_directory, LogFN);
close all;

% Display Ct(t) and fit of a single voxel
% plotSingleCTC(187, 149, 2, CTC_4D, conv_result_IRF_4D)

display('---------------------------------------------------');
display('-I- Finished Run_DCE_Perfusion execution at:');
c = clock;
disp(datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
display('---------------------------------------------------');

end

