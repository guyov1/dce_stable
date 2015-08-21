
function varargout = MainGUI(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MainGUI_OpeningFcn, ...
    'gui_OutputFcn',  @MainGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before MainGUI is made visible.
function MainGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MainGUI (see VARARGIN)

% Choose default command line output for MainGUI
handles.output = hObject;

% Update handles structure
handles.GUISaveFN=[fileparts(getComputerParams('infosfn')) filesep 'LastMainGUI.mat'];
% set(handles.DCE_Data_ListBox,'Max',3);
% set(handles.DCE_Data_ListBox,'Min',0);

% handles.destFolderBase=['/data/DCET1_' TmpName filesep];

handles.destFolderBase=[getComputerParams('basepath') 'DCEOut' filesep];

% Check if there exists a default destination folder
base_path_dir = getComputerParams('basepath');
if( exist(base_path_dir,'dir') )
    Dest_Folder_File_Path =[base_path_dir 'Default_Dest_Folder.mat'];
    if( exist(Dest_Folder_File_Path,'file') )
        default_folder = load(Dest_Folder_File_Path);
        handles.destFolderBase = default_folder.Default_dest_folder_Val;
    end
end

% If using FMRI-T9, default output folder will be D drive
% if strcmp('FMRI-T9',getenv('ComputerName'))
%     
%     %handles.destFolderBase= 'C:\Users\guyn\Desktop\DCE_OUT\';
%     handles.destFolderBase= 'D:\users\guyn\DCE_OUT\';
%     
%     % Set the Temp directory under the destination directory
%     setComputerParamM('temppath',[handles.destFolderBase 'Temp']);
% end

set(handles.Dest_Folder_TextBox,'String',handles.destFolderBase);
set(handles.list_of_paths_box,'Min',0);
set(handles.list_of_paths_box,'Max',3);

Defaults.SubSampling=1;
Defaults.nVolsToRemoveFromStart=0;
Defaults.nVolsToRemoveFromEnd=0;
Defaults.SubSecondResolution=2;
Defaults.MinFirstBolusStd=2;
Defaults.EM_Num_Of_Iterations=0;
Defaults.FMS_TolFun=1e-11;
Defaults.FMS_MaxFunEvals=10000;
Defaults.FMS_MaxIter=10000;
Defaults.MaxTDif_ForAIFSearch=3;
Defaults.MaxTDif_ForWholeVOI=6;
Defaults.Rep_MaxAroundBolus=4;
Defaults.Rep_RatioToEnd=4;
Defaults.Rep_nPerSet=2;
Defaults.Mask_thresh=0.5;
Defaults.Run_On_All=0;
Defaults.TimeDelayToMaskVeins=-0.5;
Defaults.WeightForAIFMeanVesses=0.3;
Defaults.MainCoregistration=1;
Defaults.CoregRelaxToMain=1;

% If "1" and manualArt.nii exists, take the arteries from that file, 
% take their average and make a regular calculation (we have AIF so we 
% simply use Murase to get the PK parameters) without the possibility to shift BAT.
Defaults.MakeNoBATManualArtAnalysis=0;
% If "1" and manualArt.nii exists, take the arteries from that file, calculate the parameters 
% using F Min Search on the picked arteries (instead of finding representative) and allow the possibility to shift BAT
Defaults.MakeBATManualArtAnalysis=0;
% Default mode of choosing the arteries automatically
Defaults.MakeBATAutoArtAnalysis=1;
Defaults.MakeBATPopArtAnalysis=0;
Defaults.ExtractFAs=1;
Defaults.IncludeMainInT1=0;
Defaults.UseN3OnT1=1;
Defaults.TimeMultiplier=1;
Defaults.Use_Single_M0=0;
Defaults.Calc_Gains_Diff=1;
Defaults.ThreshForRefMasks=0.99;

ParamToolTips.SubSampling=['SubSampling: Allowing to sub-sample the original data (use lower temporal resolution).' 10 'Default: 1. DO NOT CHANGE (used for high resolution data).'];
ParamToolTips.nVolsToRemoveFromEnd=['nVolsToRemoveFromEnd: Cut the last volumes of the test (sometimes the last volumes are distorted).' 10 'Default: 0'];
ParamToolTips.SubSecondResolution=['SubSecondResolution: Number of sub seconds parts for super resolution ("2" means 1/2 of a second).' 10 'Default: 2'];
ParamToolTips.MinFirstBolusStd=['MinFirstBolusStd: The minimum width of the bolus (standard deviation of the Gaussian that represents the first bolus).' 10 'Default: 2'];
ParamToolTips.EM_Num_Of_Iterations=['EM_Num_Of_Iterations: Number of iterations for the Expected Minimization algorithm which finds the optimal AIF and parameters. ).' 10 'Default: 5. (Currently not used, the algorithm uses Murase)'];
ParamToolTips.FMS_TolFun=['FMS_TolFun: Function Minimum Search''s (Matlab''s) parameter. Tolerate Function – minimal improvement for continuing the search.' 10 'Default: e^(-11)'];
ParamToolTips.FMS_MaxFunEvals=['FMS_MaxFunEvals: Number of possibilities for the F Mean Search at each step to change. Can think of it as in the case of 2-D vector f(X) ( How many 2-D  points to move to from the current one).' 10 'Default: 10000'];
ParamToolTips.FMS_MaxIter=['FMS_MaxIter: Maximal Number of iterations for FMS algorithm.' 10 'Default: 10000'];
ParamToolTips.MaxTDif_ForAIFSearch=['MaxTDif_ForAIFSearch: The possible shift in time for the AIF of the representing voxels (in seconds).' 10 'Default: 3'];
ParamToolTips.MaxTDif_ForWholeVOI=['MaxTDif_ForWholeVOI: Same as MaxTDif_ForAIFSearch, just when allowing shifting in time for all voxels in VOI (and not just representing voxels).' 10 'Default: 6'];
ParamToolTips.Rep_MaxAroundBolus=['Rep_MaxAroundBolus: Number of clusters around the bolus (for finding representing voxels).' 10 'Default: 10'];
ParamToolTips.Rep_RatioToEnd=['Rep_RatioToEnd: Number of clusters around the end of the test (for finding representing voxels).' 10 'Default: 10'];
ParamToolTips.Rep_nPerSet=['Rep_nPerSet: Number of total clusters will be MaxAroundBolus *Rep_RatioToEnd. This option will determine how many representing voxels we will choose from each cluster.' 10 'Default: 1'];
ParamToolTips.MakeNoBATManualArtAnalysis=['MakeNoBATManualArtAnalysis: If "1" and manualArt.nii exists, take the arteries from that file, take their average and make a regular calculation (we have AIF so we simply use Murase to get the PK parameters) without the possibility to shift BAT.' 10 'Default: 0'];
ParamToolTips.MakeBATManualArtAnalysis=['MakeBATManualArtAnalysis: If "1" and manualArt.nii exists, take the arteries from that file, calculate the parameters using F Min Search on the picked arteries (instead of finding representative) and allow the possibility to shift BAT.' 10 'Default: 0'];
ParamToolTips.MakeBATAutoArtAnalysis=['MakeBATAutoArtAnalysis: The default mode of choosing the arteries automatically.' 10 'Default: 1'];
ParamToolTips.ExtractFAs=['ExtractFAs: Correct the flip angles of the scan (we assume there is an error).' 10 'Default: 1'];
ParamToolTips.IncludeMainInT1=['IncludeMainInT1: Default: 1'];
ParamToolTips.UsingN3T1=['UsingN3T1: Default: 1'];
ParamToolTips.TimeMultiplier=['TimeMultiplier: Default: 1'];
ParamToolTips.Use_Single_M0=['Use_Single_M0: Enable calculating T1 using a single angle.' 10 'Default: 0'];
ParamToolTips.Calc_Gains_Diff=['Calc_Gains_Diff: Enable/disable gains calculation made by Gilad.' 10 'Default: 1'];
ParamToolTips.Mask_thresh=['Mask_thresh: Use positive values (0 to 1) to use SPM for brain masking with the given threshold.' 10 'Use negative values (0 to -1) to use FMRIB''s BET for brain masking with the given threshold (in absolute).' 10 'Default: 0.5'];
handles.ParamToolTips=ParamToolTips;

handles.Options=Defaults;
FNames=fieldnames(handles.Options);
FillTo=max(cellNumel(FNames))+2;
for i=1:numel(FNames)
    Strs{i}=[FNames{i} repmat(' ',[1 FillTo-numel(FNames{i})]) num2str(handles.Options.(FNames{i}))];
end
set(handles.Script_Params_listbox,'String',Strs);
set(handles.Script_Params_listbox,'Value',1);

tmp=getKthElement(getComputerParams('tpm'),1);
ROIsPath=[tmp(1:end-23) 'Code' filesep 'Utils' filesep 'SPM_precofigures' filesep];

checkboxAdvanced_Callback(hObject, eventdata, handles);

try
    resetGUI(hObject, handles);
catch
end

D=dir([ROIsPath '*.nii']);
ROINames=cellfun(@(x) x(1:end-4),{D.name}','UniformOutput',false);
set(handles.listboxROI,'String',['Full';ROINames]);
guidata(hObject, handles);

% UIWAIT makes MainGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = MainGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Choose_Folder_PushButton.
function Choose_Folder_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Choose_Folder_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check if there exists a default data folder
base_path_dir = getComputerParams('basepath');
if( exist(base_path_dir,'dir') )
    Data_Folder_File_Path =[base_path_dir 'Default_Data_Folder.mat'];
    if( exist(Data_Folder_File_Path,'file') )
        default_folder = load(Data_Folder_File_Path);
        handles.BaseP = default_folder.Default_data_folder_Val;
    else % Defualt case when no folder was chosen before
        handles.BaseP='/data/';
    end
end
handles.BaseP=uigetdir(handles.BaseP);

if(handles.BaseP(end) ~= filesep)
    handles.BaseP=[handles.BaseP filesep];
end

% Save the data folder as the default one for next use
base_path_dir = getComputerParams('basepath');
if( exist(base_path_dir,'dir') )
    Data_Folder_File_Path =[base_path_dir 'Default_Data_Folder.mat'];
    Default_data_folder_Val = handles.BaseP;
    save(Data_Folder_File_Path,'Default_data_folder_Val');
end


set(handles.Data_Folder_Text_Box,'String',handles.BaseP);

% JustNameFN=[handles.BaseP 'JustName.mat'];
% 
% MatFN=[handles.BaseP 'Data.mat'];
% 
% % If chosen by user, delete old Data.mat file
% if ( get(handles.Delete_Old_Data_File_CheckBox,'Value') )
%     if(exist(MatFN,'file'))
%         delete(MatFN);
%     end
% end
% 
% if(exist(MatFN,'file'))
%     handles=Read_Protocol_Action_Callback(hObject, eventdata, handles);
% end

guidata(hObject, handles);
saveGUIStatus(handles);


% --- Executes on selection change in DCE_Data_ListBox.
function DCE_Data_ListBox_Callback(hObject, eventdata, handles)
% hObject    handle to DCE_Data_ListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns DCE_Data_ListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DCE_Data_ListBox
handles.GoodSers=get(handles.DCE_Data_ListBox,'Value');
guidata(hObject, handles);
saveGUIStatus(handles);

% --- Executes during object creation, after setting all properties.
function DCE_Data_ListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DCE_Data_ListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = Process4D(hObject, eventdata, handles)

% Get indices of the data in the listbox
GoodS=get(handles.DCE_Data_ListBox,'Value');

% The following is used to get the DCE over time data
DCEMainI=getDCEMainI(handles,GoodS);

f = helpdlg([handles.ShortSeriesName ': Preprocess...']);

DCEMainDCMFN=[handles.destFolder 'DCEMain.dcm'];
copyfile(handles.ShortInfos(GoodS(DCEMainI)).Filename,DCEMainDCMFN,'f');

% Create relevant masks for the 4D data and find bolus time.
DCET1_Prepare4Df(handles.ShortInfos(GoodS(DCEMainI)),handles.destFolder,get(handles.forceP4d,'Value'),handles.Options);

% If help dialog is still open, close it
if ( exist('f','var') && ishandle(f) )
    status = close(f);
    if (~status)
        display('-W- An error occured while trying to close the help dialog!');
    end
end

% --- Executes on button press in Process_T1_PushButton.
function handles=Process_T1_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Process_T1_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function handles=ProcessT1(hObject, eventdata, handles)

% Get all the T1 data (both main and T1 calc)
GoodS=get(handles.DCE_Data_ListBox,'Value');

% Find the Main Data
DCEMainI=getDCEMainI(handles,GoodS);

% Remove the Main from the rest
GoodS2=GoodS(setdiff(1:numel(GoodS),DCEMainI));

% Get the Info of the relaxometry data
RelaxInfox=handles.ShortInfos(GoodS2(ismember({handles.ShortInfos(GoodS2).Class},{'DCE','T1MAP'})));

% Check if N3 was chosen
DoN3 = get(handles.N3_B1_CheckBox,'Value');

% ASK GILAD - what is the purpose of DoGlobal and DoDeviations?
% ANSWER -
% DoDeviations ï¿½ Correction for flip angles (exp. ï¿½ machine says its 20 degrees while the real one is 20.2).
%                                     Implemented under CalcRelaxForVolNFA.m. Currently not used.
%
% DoGlobal ï¿½ N3 fixes B1 non-uniformity by a factor for each voxel.
%                          Gilad picks a point in white matter to be the reference point (implemented in DCET1_RelaxForSubjectf.m).
DoGlobal     = get(handles.Global_FA_Factor_CheckBox,'Value');
DoDeviations = get(handles.FA_Deviations_CheckBox,'Value');

% See if we chose to force calc T1 (if not, calculation is skipped if was previously done)
CalcForce=get(handles.forcet1,'Value');

% Get from GUI to what to co-register T1 maps
CrgToS=get(handles.t1coregto,'String');
CrgTo=CrgToS{get(handles.t1coregto,'Value')};

% Find the main data file
GoodS=get(handles.DCE_Data_ListBox,'Value');
DCEMainI=getDCEMainI(handles,GoodS);
MainInfo=handles.ShortInfos(GoodS(DCEMainI));

% Starting T1 calculation
f = helpdlg([handles.ShortSeriesName ': T1 calculation... ']);

% Calculate T1 maps
DCET1_RelaxForSubjectf(RelaxInfox,handles.destFolder,DoN3,DoGlobal,DoDeviations,CalcForce,CrgTo,handles.Options,MainInfo);

% If help dialog is still open, close it
if ( exist('f','var') && ishandle(f) )
    
    status = close(f);
    if (~status)
        display('-W- An error occured while trying to close the help dialog!');
    end
end

% --- Executes on button press in Process_T2_PushButton.
function Process_T2_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Process_T2_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Create_CTCs_PushButton.
function handles=Create_CTCs_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Create_CTCs_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function handles=CreateCTCs(hObject, eventdata, handles)

% Check if N3 was chosen
DoN3 = get(handles.N3_B1_CheckBox,'Value');

% DoDeviations ï¿½ Correction for flip angles (exp. ï¿½ machine says its 20 degrees while the real one is 20.2).
%                                     Implemented under CalcRelaxForVolNFA.m. Currently not used.
% DoGlobal ï¿½ N3 fixes B1 non-uniformity by a factor for each voxel.
%                          Gilad picks a point in white matter to be the reference point (implemented in DCET1_RelaxForSubjectf.m).
DoGlobal           = get(handles.Global_FA_Factor_CheckBox,'Value');
DoDeviations       = get(handles.FA_Deviations_CheckBox,'Value');

% See if user chose to forcly calculate CTC
CalcForce=get(handles.forcectc,'Value');

% Get all the T1 data (both main and T1 calc)
GoodS=get(handles.DCE_Data_ListBox,'Value');

% Find the main data file
DCEMainI=getDCEMainI(handles,GoodS);

% Remove the Main from the rest
GoodS2=GoodS(setdiff(1:numel(GoodS),DCEMainI));
Additional_T1_Maps_Time_Stamps = {handles.ShortInfos(GoodS2(ismember({handles.ShortInfos(GoodS2).Class},{'DCE','T1MAP'}))).SeriesDateTime};

% UnderSampling=str2num(get(handles.Under_Sampling_EditBox,'String'));
% See if user chose UnderSampling parameter -
% Allowing to sub-sample  the original data and then use super resolution
UnderSampling=handles.Options.SubSampling;

% Get from GUI to what to co-register T1 maps
CrgToS=get(handles.t1coregto,'String');
CrgTo=CrgToS{get(handles.t1coregto,'Value')};

% Calculate CTC
f = helpdlg([handles.ShortSeriesName ': CTCs calculation...']);
% Check for additional DESPOT1 maps for 4D data
DCET1_CTCf(handles.ShortInfos(GoodS(DCEMainI)),...
           Additional_T1_Maps_Time_Stamps,handles.destFolder,DoN3,DoGlobal,DoDeviations,CalcForce,CrgTo,handles.Options);

% If help dialog is still open, close it
if ( exist('f','var') && ishandle(f) )
    
    status = close(f);
    if (~status)
        display('-W- An error occured while trying to close the help dialog!');
    end
end

% --- Executes on button press in DCE_Analysis_PushButton.
function DCE_Analysis_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to DCE_Analysis_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function ArtFind(hObject, eventdata, handles)
WorkingP=handles.destFolder;
% UnderSampling=str2num(get(handles.Under_Sampling_EditBox,'String'));
UnderSampling=handles.Options.SubSampling;
if(UnderSampling==1)
    USStr='';
else
    USStr=['_' num2str(UnderSampling)];
end
CTCFN=[WorkingP 'AfterCTC' USStr '.mat'];
disp('Loading..');
load(CTCFN);

Options=handles.Options;

GoodS=get(handles.DCE_Data_ListBox,'Value');
DCEMainI=getDCEMainI(handles,GoodS);

CurMainDCEInfo=handles.ShortInfos(GoodS(DCEMainI));
f = helpdlg([handles.ShortSeriesName ': Finding represntative voxels...']);
disp('DCET1_PKf..');
DCET1_PKf;
disp('AIFFinderTry..');
ROIStrs=get(handles.listboxROI,'String');
ROIStr=ROIStrs{get(handles.listboxROI,'value')};
ArtFinder;
% If help dialog is still open, close itit
if ( exist('f','var') && ishandle(f) )
    status = close(f);
    if (~status)
        display('-W- An error occured while trying to close the help dialog!');
    end
end

function AIFFind(hObject, eventdata, handles)
WorkingP=handles.destFolder;
% UnderSampling=str2num(get(handles.Under_Sampling_EditBox,'String'));
UnderSampling=handles.Options.SubSampling;
if(UnderSampling==1)
    USStr='';
else
    USStr=['_' num2str(UnderSampling)];
end
CTCFN=[WorkingP 'AfterCTC' USStr '.mat'];
disp('Loading..');
load(CTCFN);

Options=handles.Options;

GoodS=get(handles.DCE_Data_ListBox,'Value');
DCEMainI=getDCEMainI(handles,GoodS);

CurMainDCEInfo=handles.ShortInfos(GoodS(DCEMainI));
f = helpdlg([handles.ShortSeriesName ': AIF extraction...']);
disp('DCET1_PKf..');
DCET1_PKf;
disp('AIFFinderTry..');
AIFFinder;
% If help dialog is still open, close it
if ( exist('f','var') && ishandle(f) )
    status = close(f);
    if (~status)
        display('-W- An error occured while trying to close the help dialog!');
    end
end

function WholeROICompute(hObject, eventdata, handles)
WorkingP=handles.destFolder;
% UnderSampling=str2num(get(handles.Under_Sampling_EditBox,'String'));
UnderSampling=handles.Options.SubSampling;
if(UnderSampling==1)
    USStr='';
else
    USStr=['_' num2str(UnderSampling)];
end
CTCFN=[WorkingP 'AfterCTC' USStr '.mat'];
disp('Loading..');
load(CTCFN);

Options=handles.Options;

GoodS=get(handles.DCE_Data_ListBox,'Value');
DCEMainI=getDCEMainI(handles,GoodS);

CurMainDCEInfo=handles.ShortInfos(GoodS(DCEMainI));

f = helpdlg([handles.ShortSeriesName ': Whole ROI computation...']);
disp('DCET1_PKf..');
% Get clustering, bolus time and noise...
DCET1_PKf;
disp('WholeROICompute..');

WorkingP=handles.destFolder;

% Calculate parameters for the entire brain (all voxels)
DCET1_WholeVolCompute;

% If help dialog is still open, close it
if ( exist('f','var') && ishandle(f) )
    
    status = close(f);
    if (~status)
        display('-W- An error occured while trying to close the help dialog!');
    end
end

% --- Executes on button press in N3_B1_CheckBox.
function N3_B1_CheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to N3_B1_CheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of N3_B1_CheckBox


% --- Executes on button press in Global_FA_Factor_CheckBox.
function Global_FA_Factor_CheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to Global_FA_Factor_CheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Global_FA_Factor_CheckBox


% --- Executes on button press in FA_Deviations_CheckBox.
function FA_Deviations_CheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to FA_Deviations_CheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FA_Deviations_CheckBox


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in Choose_Dest_Folder_PushButton.
function Choose_Dest_Folder_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Choose_Dest_Folder_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(exist(handles.destFolderBase,'dir'))
    TmpP=uigetdir(handles.destFolderBase);
    % Set the Temp directory under the destination directory
    setComputerParamM('temppath',[TmpP 'Temp']);
    
else
    TmpP=uigetdir('/data/');
end
if(numel(TmpP)<4)
    return;
end
if(TmpP(end)~=filesep)
    TmpP=[TmpP filesep];
end

handles.destFolderBase=TmpP;

% Save the destination as the default one for next use
base_path_dir = getComputerParams('basepath');
if( exist(base_path_dir,'dir') )
    Dest_Folder_File_Path =[base_path_dir 'Default_Dest_Folder.mat'];
    Default_dest_folder_Val = handles.destFolderBase;
    % Remove old default file
    %if(exist(Dest_Folder_File_Path,'file'))
    %    delete(Dest_Folder_File_Path);
    %end
    save(Dest_Folder_File_Path,'Default_dest_folder_Val');
end
   

set(handles.Dest_Folder_TextBox,'String',handles.destFolderBase);
guidata(hObject, handles);

% --- Executes on button press in forceP4d.
function forceP4d_Callback(hObject, eventdata, handles)
% hObject    handle to forceP4d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of forceP4d


% --- Executes on button press in SliceIntensity_BolusTime_PushButton.
function SliceIntensity_BolusTime_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to SliceIntensity_BolusTime_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
system(['gthumb ' handles.destFolder 'SlicesIntensityAndBolusTime.jpg']);

% --- Executes on button press in Coreged_Mid_Slice_PushButton.
function Coreged_Mid_Slice_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Coreged_Mid_Slice_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
system(['gthumb ' handles.destFolder 'CoregedMidSlice.jpg']);


% --- Executes on button press in MeanVolume_PushButton.
function MeanVolume_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to MeanVolume_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mricronx([handles.destFolder 'DCEMean.nii']);


% --- Executes on selection change in t1coregto.
function t1coregto_Callback(hObject, eventdata, handles)
% hObject    handle to t1coregto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns t1coregto contents as cell array
%        contents{get(hObject,'Value')} returns selected item from t1coregto


% --- Executes during object creation, after setting all properties.
function t1coregto_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t1coregto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in forcet1.
function forcet1_Callback(hObject, eventdata, handles)
% hObject    handle to forcet1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of forcet1

% --- Executes on button press in Do_Debug_CheckBox.
function Do_Debug_CheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to Do_Debug_CheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Do_Debug_CheckBox


% --- Executes on button press in forcectc.
function forcectc_Callback(hObject, eventdata, handles)
% hObject    handle to forcectc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of forcectc

function saveGUIStatus(handles)
% Method for saving GUI elements properties
% obj  the object that sends the data
% handles  the GUI handles structure
handlesFields = fieldnames(handles); % extract all the field names in handles
tmpdata = cell(size(handlesFields)); % make a tmp data cell
GUIprop = cell2struct(tmpdata, handlesFields, 1); % make a structure with the same field names as handles
for(i = 1:length(handlesFields)) % index all the fields
    test = ishandle(handles.(handlesFields{i})); % test if the name of field is a handle reference
    if(test)
        GUIprop.(handlesFields{i}) = get(findobj('Tag', handlesFields{i})); % save the properties
    else
        GUIprop.(handlesFields{i}) = handles.(handlesFields{i}); %rmfield(obj.GUIprop, handlesFields{i}); % remove the element because is not a GUI handle
    end
end
save(handles.GUISaveFN,'GUIprop');

%  end code here
% The application is using object oriented programing but one can disregard the obj. structure like reference. The code above basically saves only the widgets from the handles structure. By using the get function with only the tag name of the element I extract all the properties of the widget.
% Now I will load that data saved previously in the GUIprop structure. The code here:
% - start code here

function handles=resetGUI(hObject, handles)
% Method for reseting GUI elements properties
% obj  the object that sends the data
% handles  the GUI handles structure
GUIprop=load(handles.GUISaveFN);
GUIprop=GUIprop.GUIprop;

handlesFields = fieldnames(GUIprop); % extract all the field names in handles
for(i = 2:length(handlesFields)) % index all the fields (do not take figure1 properties)
    test = isfield(handles,handlesFields{i}) && any(ishandle(handles.(handlesFields{i}))); % test if the name of field is a handle reference
    if(test)
        setableProps = set(findobj('Tag', handlesFields{i})); % extract the user accessible properties
        if(isempty(setableProps))
            continue;
        end
        setableNames = fieldnames(setableProps); % extract only the names of properties
        allProps = fieldnames(GUIprop.(handlesFields{i})); % all the properties of the element
        [tf, loc] = ismember(allProps, setableNames); % map the properties to update
        valuesOriginal = struct2cell(GUIprop.(handlesFields{i})); % extract the values
        values = valuesOriginal(tf); % filter the accessible values
        for j=1:numel(setableNames)
            if(~ismember(setableNames{j},{'Parent','Position','Visible'}))
                set(findobj('Tag', handlesFields{i}), setableNames{j}, values{j}); % save the properties
            end
        end
    else
        handles.(handlesFields{i})=GUIprop.(handlesFields{i});
    end
end
% guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Choose_Folder_PushButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Choose_Folder_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in Read_Protocol_Action.
function handles=Read_Protocol_Action_Callback(hObject, eventdata, handles)
% hObject    handle to Read_Protocol_Action (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(~isfield(handles,'BaseP'))
    return;
end
CurBaseP=handles.BaseP;
CurBaseP=strrep(CurBaseP,[filesep filesep],filesep);
F=find(CurBaseP=='|');
if(~isempty(F))
    CurBaseP=CurBaseP(F(end)+2:end);
end

if (filesep ~='/') % If not Unix
    if(CurBaseP(1)==filesep  && CurBaseP(2) ~=filesep)
        CurBaseP=[filesep CurBaseP];
    end
end

JustNameFN=[CurBaseP 'JustName.mat'];
if(exist(JustNameFN,'file'))
    load(JustNameFN);
else
    ShortInfos = ReadNewScans(CurBaseP);
    ShortSeriesName=[ToShortName(ShortInfos(end).Name) '_' ShortInfos(end).SeriesDate];
    save(JustNameFN,'ShortSeriesName');
end
% Remove any spaces
ShortSeriesName = regexprep(ShortSeriesName,'\W','');

handles.ShortSeriesName=ShortSeriesName;
handles.destFolder=[handles.destFolderBase handles.ShortSeriesName filesep];
[s, mess, messid] =mkdir(handles.destFolder);

ShortInfosFN=[handles.destFolder 'ShortInfos.mat'];
if(exist('ShortInfos','var'))
    save(ShortInfosFN,'ShortInfos');
else
    if(exist(ShortInfosFN,'file'))
        load(ShortInfosFN);
    else
        ShortInfos = ReadNewScans(CurBaseP);
        save(ShortInfosFN,'ShortInfos');
    end
end
handles.ShortInfos=ShortInfos;

set(handles.Data_Folder_Text_Box,'String',handles.ShortSeriesName);
Descriptions=strcat({handles.ShortInfos.SeriesDescription},'_');
Descriptions=strcat(Descriptions,{handles.ShortInfos.Class});
Descriptions=strrep(Descriptions,'_',' ');

GoodS=intersect(get(handles.DCE_Data_ListBox,'Value'),1:numel(handles.ShortInfos));
Tmp=[handles.ShortInfos(GoodS).ImagesInAcquisition];
Tmp(Tmp==25)=2500; % For PHILIPS now
DCEMainIs=find(Tmp>50); % Siemens

for i=DCEMainIs
    Descriptions{i}=[Descriptions{i} ' ' num2str(handles.ShortInfos(i).FlipAngle) ' ' num2str(handles.ShortInfos(i).ImagesInAcquisition)];
end

set(handles.DCE_Data_ListBox,'String',Descriptions);
set(handles.DCE_Data_ListBox,'Max',3);
set(handles.DCE_Data_ListBox,'Min',0);

set(handles.Dest_Folder_TextBox,'String',handles.destFolder);

handles.GoodSers=find(ismember({handles.ShortInfos.Class},{'T1MAP','DCE'})); % 'FIESTA'

% GUYN - DATA CHANGE - commented to select all T1 data and not just up to Main data

% Select just the data we need
set(handles.DCE_Data_ListBox,'Value',handles.GoodSers);

guidata(hObject, handles);
saveGUIStatus(handles);


% --- Executes on button press in Do_All_Actions.
function handles=Do_All_Actions_Callback(hObject, eventdata, handles)
% hObject    handle to Do_All_Actions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=Read_Protocol_Action_Callback(hObject, eventdata, handles);
% Get string of all possible stages
Strs=get(handles.Stages_List_Box,'String');



% If user chose run on all, do all stages. Else, do selected ones.
if (handles.Options.Run_On_All)
    ToDo=Strs;
else % run on selected
    V=get(handles.Stages_List_Box,'Value');
    ToDo=Strs(V);
end

WorkingP=handles.destFolder;


% Check if user chose debug mode
DoDebug=get(handles.Do_Debug_CheckBox,'Value');

% In DEBUG mode, don't use try-catch
if(DoDebug)
    if(~exist(handles.destFolderBase,'dir'))
        mkdir(handles.destFolderBase);
    end
    if(ismember('Process 4D',ToDo))
        handles=Process4D(hObject, eventdata, handles);
        MakeReport; % Create PDF Report
    end
    if(ismember('Process T1',ToDo))
        handles=ProcessT1(hObject, eventdata, handles);
        MakeReport; % Create PDF Report
    end
    if(ismember('Create CTCs',ToDo))
        handles=CreateCTCs(hObject, eventdata, handles);
        MakeReport; % Create PDF Report
    end
    if(ismember('Find arteries',ToDo))
        ArtFind(hObject, eventdata, handles);
        MakeReport; % Create PDF Report
    end
    if(ismember('Find AIF',ToDo))
        AIFFind(hObject, eventdata, handles);
        MakeReport; % Create PDF Report
    end
    if(ismember('Whole ROI computation',ToDo))
        WholeROICompute(hObject, eventdata, handles);
        MakeReport; % Create PDF Report
    end
else % In regular mode, use try-catch
    try
        if(~exist(handles.destFolderBase,'dir'))
            mkdir(handles.destFolderBase);
        end
        if(ismember('Process 4D',ToDo))
            handles=Process4D(hObject, eventdata, handles);
            MakeReport; % Create PDF Report
        end
        if(ismember('Process T1',ToDo))
            handles=ProcessT1(hObject, eventdata, handles);
            MakeReport; % Create PDF Report
        end
        if(ismember('Create CTCs',ToDo))
            handles=CreateCTCs(hObject, eventdata, handles);
            MakeReport; % Create PDF Report
        end
        if(ismember('Find arteries',ToDo))
            ArtFind(hObject, eventdata, handles);
            MakeReport; % Create PDF Report
        end
        if(ismember('Find AIF',ToDo))
            AIFFind(hObject, eventdata, handles);
            MakeReport; % Create PDF Report
        end
        if(ismember('Whole ROI computation',ToDo))
            WholeROICompute(hObject, eventdata, handles);
            MakeReport; % Create PDF Report
        end
    % Catch error struct when needed    
    catch ME
        
        % Add to log file the error details
        StackSize=numel(ME.stack);
        AddToLog(WorkingP,'zzzz_0a','\\section*{Error:}',[],3);
        AddToLog(WorkingP,'zzzz_0b',CorrectStrForLog(ME.identifier));
        AddToLog(WorkingP,'zzzz_0c',CorrectStrForLog(ME.message));
        for ii=1:StackSize
            AddToLog(WorkingP,['zzzz_' num2str(ii) 'a'],CorrectStrForLog(ME.stack(ii).file));
            AddToLog(WorkingP,['zzzz_' num2str(ii) 'b'],CorrectStrForLog([ME.stack(ii).name ' ' num2str(ME.stack(ii).line)]));
        end
        
        MakeReport; % Create PDF Report
    end
end

% --- Executes on selection change in list_of_paths_box.
function list_of_paths_box_Callback(hObject, eventdata, handles)
% hObject    handle to list_of_paths_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_of_paths_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_of_paths_box
if(~isfield(handles,'LastTimeClickedList'))
    handles.LastTimeClickedList=now;
    guidata(hObject, handles);
    return;
end
TD=now-handles.LastTimeClickedList;
handles.LastTimeClickedList=now;
if(TD<3e-6)
    S=get(handles.list_of_paths_box,'String');
    handles.BaseP=S{get(handles.list_of_paths_box,'Value')};
    
    set(handles.Data_Folder_Text_Box,'String',handles.BaseP);
    
    %     MatFN=[handles.BaseP 'Data.mat'];
    %     if(exist(MatFN,'file'))
    
    handles=Read_Protocol_Action_Callback(hObject, eventdata, handles);
    %     end
end
guidata(hObject, handles);
saveGUIStatus(handles);

% --- Executes during object creation, after setting all properties.
function list_of_paths_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_of_paths_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Add_To_List_PushButton.
function Add_To_List_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Add_To_List_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurPs=get(handles.list_of_paths_box,'String');
if(~iscell(CurPs))
    CurPs={CurPs};
    CurPs=CurPs(cellNumel(CurPs)>0);
end
CurPs=[CurPs; handles.BaseP];
set(handles.list_of_paths_box,'String',CurPs);
handles.CurPs=CurPs;
guidata(hObject, handles);
saveGUIStatus(handles);

% --- Executes on button press in Run_On_All_PushButton.
function Run_On_All_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Run_On_All_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

s=warning('off','MATLAB:MKDIR:DirectoryExists');

handles.ListMode=true;

Strs=get(handles.list_of_paths_box,'String');

if(~iscell(Strs))
    Strs={Strs};
    Strs=Strs(cellNumel(Strs)>0);
end

ErrorList=[];
MEs=cell(1);
ProbPs=cell(1);
CurPs=Strs;

for i=1:numel(Strs)
    F=find(Strs{i}=='|');
    if(~isempty(F))
        CurPs{i}=Strs{i}(F(end)+2:end);
    end
end

% Get debug button selection from user
DoDebug=get(handles.Do_Debug_CheckBox,'Value');
% Mark "Run On All" flag
handles.Options.Run_On_All=1;

for i=1:numel(CurPs)
    if(DoDebug)
        handles.BaseP=CurPs{i};
        handles=Do_All_Actions_Callback(hObject, eventdata, handles);
    else
        try
            handles.BaseP=CurPs{i};
            handles=Do_All_Actions_Callback(hObject, eventdata, handles);
        catch ME
            disp(['Failed ' handles.BaseP]);
            ErrorList=[ErrorList i];
            MEs{i}=ME;
            ProbPs{i}=handles.BaseP;
        end
    end
end
handles.ListMode=false;
guidata(hObject, handles);
saveGUIStatus(handles);
warning(s);
if(isempty(ErrorList))
    disp('No errors found.');
else
    disp('Errors:');
    disp(ErrorList);
end
assignin('base','MEs',MEs);
assignin('base','ProbPs',ProbPs);
assignin('base','ErrorList',ErrorList);
disp('Finished List');


% --- Executes on button press in Remove_From_List_PushButton.
function Remove_From_List_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Remove_From_List_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Strs=get(handles.list_of_paths_box,'String');
V=get(handles.list_of_paths_box,'Value');
if(~iscell(Strs))
    Strs={Strs};
    Strs=Strs(cellNumel(Strs)>0);
end
Others=setdiff(1:numel(Strs),V);
set(handles.list_of_paths_box,'String',Strs(Others));
set(handles.list_of_paths_box,'Value',[]);

% --- Executes on button press in Clear_List_PushButton.
function Clear_List_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Clear_List_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.list_of_paths_box,'String',cell(0));
handles.CurPs=cell(0);
guidata(hObject, handles);
saveGUIStatus(handles);

% --- Executes on button press in Add_Base_Path_PushButton.
function Add_Base_Path_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Add_Base_Path_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DD=dirdfs(handles.BaseP);
DDC={DD.name};
DDC2=DDC(strhas(DDC,'\.dcm'));
for i=1:numel(DDC2)
    S=find(DDC2{i}==filesep);
    if(numel(S)>1)
        DDC3{i}=DDC2{i}(1:S(end-1));
    else
        DDC3{i}=DDC2{i}(1:S(end));
    end
end
UDDC3=unique(DDC3)';
FUDDC3=MultiConcat([handles.BaseP],UDDC3)';
CurPs=get(handles.list_of_paths_box,'String');
if(~iscell(CurPs))
    CurPs={CurPs};
    CurPs=CurPs(cellNumel(CurPs)>0);
end
CurPs=[CurPs; FUDDC3];
set(handles.list_of_paths_box,'String',CurPs);
handles.CurPs=CurPs;
guidata(hObject, handles);
saveGUIStatus(handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over list_of_paths_box.
function list_of_paths_box_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to list_of_paths_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=6;


% --- Executes on key press with focus on list_of_paths_box and none of its controls.
function list_of_paths_box_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to list_of_paths_box (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
a=7;


% --- Executes on button press in Explore_Results_PushButton.
function Explore_Results_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Explore_Results_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=Read_Protocol_Action_Callback(hObject, eventdata, handles);
if(~isfield(handles,'destFolder'))
    msgbox('Choose (double-click) a scan');
    return;
end

WorkingP=handles.destFolder;
% UnderSampling=str2num(get(handles.Under_Sampling_EditBox,'String'));
UnderSampling=handles.Options.SubSampling;
if(UnderSampling==1)
    USStr='';
else
    USStr=['_' num2str(UnderSampling)];
end
CTCFN=[WorkingP 'AfterCTC' USStr '.mat'];
disp('Loading..');
load(CTCFN);

GoodS=get(handles.DCE_Data_ListBox,'Value');
DCEMainI=getDCEMainI(handles,GoodS);

CurMainDCEInfo=handles.ShortInfos(GoodS(DCEMainI));
disp('DCET1_PKf..');
WorkingP=handles.destFolder;
PrepareFN=[WorkingP 'AfterPrepare4D.mat'];
CTCFN=[WorkingP 'AfterCTC' USStr '.mat'];
MeanFN=[WorkingP 'DCEMean.nii'];
DCET1_PKf;
disp('ExploreDCEResults..');
Options=handles.Options;
ExploreDCEResults;


% --- Executes on button press in Run_After_AIF_On_All_PushButton.
function Run_After_AIF_On_All_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Run_After_AIF_On_All_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Run_On_Selected_Button.
function Run_On_Selected_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Run_On_Selected_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Run on selected

% Disable "Directory Exists" warnings
s=warning('off','MATLAB:MKDIR:DirectoryExists');

handles.ListMode=true;
Strs=get(handles.list_of_paths_box,'String');

if(~iscell(Strs))
    Strs={Strs};
    Strs=Strs(cellNumel(Strs)>0);
else
    Strs=Strs(get(handles.list_of_paths_box,'Value'));
end

ErrorList=[];
MEs=cell(1);
ProbPs=cell(1);
CurPs=Strs;
for i=1:numel(Strs)
    F=find(Strs{i}=='|');
    if(~isempty(F))
        CurPs{i}=Strs{i}(F(end)+2:end);
    end
end

% Indication whether to use catch and try or not.
DoDebug=get(handles.Do_Debug_CheckBox,'Value');
% Unmark "Run On All" flag, becuase user chose run on selected
handles.Options.Run_On_All=0;

for i=1:numel(CurPs)
    
    if(DoDebug)
        
        handles.BaseP=CurPs{i};
        % The actual functioning
        handles=Do_All_Actions_Callback(hObject, eventdata, handles);
        
    else
        
        try
            handles.BaseP=CurPs{i};
            % The actual functioning
            handles=Do_All_Actions_Callback(hObject, eventdata, handles);
        catch ME
            disp(['Failed ' handles.BaseP]);
            ErrorList=[ErrorList i];
            MEs{i}=ME;
            ProbPs{i}=handles.BaseP;
        end
        
    end
end

handles.ListMode=false;
guidata(hObject, handles);
saveGUIStatus(handles);
warning(s);

if(isempty(ErrorList))
    disp('No errors found.');
else
    disp('Errors:');
    disp(ErrorList);
end
assignin('base','MEs',MEs);
assignin('base','ProbPs',ProbPs);
assignin('base','ErrorList',ErrorList);
disp('Finished selected');

% --- Executes on selection change in Stages_List_Box.
function Stages_List_Box_Callback(hObject, eventdata, handles)
% hObject    handle to Stages_List_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Stages_List_Box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Stages_List_Box

% --- Executes during object creation, after setting all properties.
function Stages_List_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Stages_List_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in List_From_DataBase_PushButton.
function List_From_DataBase_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to List_From_DataBase_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = helpdlg('Reading from database...');

DCE_Infos = getSpecificInfos('Class','DCE');
%Images_num = [DCE_Infos.ImagesInAcquisition];
%Needed_idx = find(Images_num(1,:)>100);
%BB=DCE_Infos(Needed_idx);
DCE_Main_Info = DCE_Infos([DCE_Infos.ImagesInAcquisition]>100);

for i=1:numel(DCE_Main_Info)
    
    Tmp = find(DCE_Main_Info(i).Path==filesep);
    Before{i} = DCE_Main_Info(i).Path(1:Tmp(end));
    
end

CurPs=unique(Before)';
set(handles.list_of_paths_box,'String',CurPs);
handles.CurPs=CurPs;

% If help dialog is still open, close it
if ( exist('f','var') && ishandle(f) )
    
    status = close(f);
    if (~status)
        display('-W- An error occured while trying to close the help dialog!');
    end
end

f = helpdlg('Checking scans status...');

LastStep=zeros(1,numel(CurPs));
LastDate=zeros(1,numel(CurPs));
StepsTitles={'Read protocol','Preprocessed','T1 calculated','CTC extracted','AIF found','3d maps calculated'};

for i=1:numel(CurPs)
    %     disp(['Checking status for ' num2str(i)]);
    
    MatFN=[handles.CurPs{i} 'Data.mat'];
    
    if(exist(MatFN,'file'))
        load(MatFN);
    else
        ShortInfos=getSpecificInfos('xPath',handles.CurPs{i});
        if(isempty(ShortInfos))
            ShortInfos = ReadNewScans(handles.CurPs{i});
        end
        save(MatFN,'ShortInfos');
    end
    
    DCEMainI=getDCEMainI(handles,GoodS);
    
    if(isempty(DCEMainI))
        continue;
    end
    
    CurShortSeriesName=[ToShortName(ShortInfos(DCEMainI(end)).Name) '_' ShortInfos(DCEMainI(end)).SeriesDate];
    CurdestFolder=[handles.destFolderBase CurShortSeriesName filesep];
    
    StepsFNs{1}=[CurdestFolder 'AfterPrepare4D.mat'];
    StepsFNs{2}=[CurdestFolder 'Relaxometry' filesep 'NFARes.mat'];
    StepsFNs{3}=[CurdestFolder 'AfterCTC.mat'];
    StepsFNs{4}=[CurdestFolder 'AIFFindData.mat'];
    StepsFNs{5}=[CurdestFolder 'PKM3D.mat'];
    
    LastStep(i)=1;
    LastDate(i)=now;
    for j=1:5
        if(exist(StepsFNs{j},'file'))
            LastStep(i)=j+1;
            D=dir(StepsFNs{j});
            LastDate(i)=D.datenum;
        end
    end
    HasManualArtStr='';
    if(exist([CurdestFolder 'manualArt.nii'],'file'))
        HasManualArtStr='has manualArt ';
    end
    
    %     Strs{i}=[num2str(LastStep(i)) ' | ' datestr(LastDate(i)) ' | ' CurPs{i}];
    Strs{i}=[HasManualArtStr StepsTitles{LastStep(i)} ' | ' datestr(LastDate(i)) ' | ' CurPs{i}];
end

set(handles.list_of_paths_box,'String',Strs);
% If help dialog is still open, close it
if ( exist('f','var') && ishandle(f) )
    
    status = close(f);
    if (~status)
        display('-W- An error occured while trying to close the help dialog!');
    end
end
guidata(hObject, handles);
saveGUIStatus(handles);


% --- Executes on button press in Force_Read_Prot_CheckBox.
function Force_Read_Prot_CheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to Force_Read_Prot_CheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Force_Read_Prot_CheckBox


% --- Executes on button press in Add_To_DataBase_PushButton.
function Add_To_DataBase_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Add_To_DataBase_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddNewScans(handles.BaseP,1);
List_From_DataBase_PushButton_Callback(hObject, eventdata, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
DestFolderBase=handles.destFolderBase;
try
    saveGUIStatus(handles);
%     save(handles.GUISaveFN,'DestFolderBase');
catch
end
delete(hObject);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% DestFolderBase=handles.destFolderBase;
% save(handles.GUISaveFN,'DestFolderBase');
% delete(hObject);



function Under_Sampling_EditBox_Callback(hObject, eventdata, handles)
% hObject    handle to Under_Sampling_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Under_Sampling_EditBox as text
%        str2double(get(hObject,'String')) returns contents of Under_Sampling_EditBox as a double
CurV=get(handles.Script_Params_listbox,'Value');
Tmp=get(handles.Script_Params_listbox,'String');
for i=1:numel(Tmp)
    FNames{i}=getKthElement(regexp(Tmp{i},'\W+','split'),1);
end
% FNames=fieldnames(handles.Options);
Str=get(handles.Under_Sampling_EditBox,'String');
Val=str2double(Str);
if(isnan(Val))
    msgbox('Please enter numeric value only.');
    return;
end
i=CurV;
handles.Options.(FNames{i})=Val;
FillTo=max(cellNumel(FNames))+2;
Strs=get(handles.Script_Params_listbox,'String');
Strs{i}=[FNames{i} repmat(' ',[1 FillTo-numel(FNames{i})]) num2str(handles.Options.(FNames{i}))];
set(handles.Script_Params_listbox,'String',Strs);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Under_Sampling_EditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Under_Sampling_EditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Script_Params_listbox.
function Script_Params_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to Script_Params_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Script_Params_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Script_Params_listbox
CurV=get(handles.Script_Params_listbox,'Value');
Tmp=get(handles.Script_Params_listbox,'String');
for i=1:numel(Tmp)
    FNames{i}=getKthElement(regexp(Tmp{i},'\W+','split'),1);
end
% FNames=fieldnames(handles.Options);
set(handles.Under_Sampling_EditBox,'String',num2str(handles.Options.(FNames{CurV})));
if(numel(CurV)>0 && isfield(handles.ParamToolTips,FNames{min(CurV)}))
    set(handles.Script_Params_listbox,'Tooltip',handles.ParamToolTips.(FNames{min(CurV)}));
end

% --- Executes during object creation, after setting all properties.
function Script_Params_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Script_Params_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Read_Protocol_PushButton.
function Read_Protocol_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Read_Protocol_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Do_All_Actions_PushButton.
function Do_All_Actions_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Do_All_Actions_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in t1coregto.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to t1coregto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns t1coregto contents as cell array
%        contents{get(hObject,'Value')} returns selected item from t1coregto


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t1coregto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Delete_Old_Data_File_CheckBox.
function Delete_Old_Data_File_CheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to Delete_Old_Data_File_CheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Delete_Old_Data_File_CheckBox


% --- Executes on button press in Forget_Old_Info_PushButton.
% Remove the Database directory (which contains the info) and run DCEInit again
function Forget_Old_Info_PushButton_Callback(hObject, eventdata, handles)

% Get Matlab user name
TmpName= license('inuse');
TmpName=TmpName(1).user;

% Remove old CompParams.mat
if (exist('CompParams.mat', 'file'))
    delete( which('CompParams.mat') );
end

% Set base path
BaseBaseP=[pwd filesep];


% Delete old user folder (create again in DCEInit)
if ( exist([BaseBaseP TmpName], 'file') )
    outputString = sprintf('-I- Removing Directory: %s ', [BaseBaseP TmpName]);
    display(outputString);
    [SUCCESS,MESSAGE,MESSAGEID] = rmdir([BaseBaseP TmpName],'s');
    if (~SUCCESS)
        display(['-W- Could not remove:' BaseBaseP TmpName]);
        display(MESSAGE);
    end
    
end

% After removing run DCEInit again
display('-I- Running DCEInit again!');
DCEInit;


% hObject    handle to Forget_Old_Info_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in OpenData.
function OpenData_Callback(hObject, eventdata, handles)
% hObject    handle to OpenData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~isfield(handles,'BaseP'))
    return;
end
CurBaseP=handles.BaseP;
CurBaseP=strrep(CurBaseP,[filesep filesep],filesep);
F=find(CurBaseP=='|');
if(~isempty(F))
    CurBaseP=CurBaseP(F(end)+2:end);
end

if (filesep ~='/') % If not Unix
    if(CurBaseP(1)==filesep  && CurBaseP(2) ~=filesep)
        CurBaseP=[filesep CurBaseP];
    end
end
Explore(CurBaseP);

% --- Executes on button press in OpenDest.
function OpenDest_Callback(hObject, eventdata, handles)
% hObject    handle to OpenDest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Explore(handles.destFolder);


% --- Executes on button press in pushManArt.
function pushManArt_Callback(hObject, eventdata, handles)
% hObject    handle to pushManArt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RepVoxGUI(handles.destFolder);


% --- Executes on button press in pushManAIF.
function pushManAIF_Callback(hObject, eventdata, handles)
% hObject    handle to pushManAIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AIFGUI(handles.destFolder);


% --- Executes on button press in checkboxAdvanced.
function checkboxAdvanced_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAdvanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxAdvanced
AdvancedOptions={'SubSampling','EM_Num_Of_Iterations','FMS_TolFun','FMS_MaxFunEvals','FMS_MaxIter','MakeBATManualArtAnalysis','MakeBATAutoArtAnalysis','Calc_Gains_Diff'};

if(get(handles.checkboxAdvanced,'Value'))
    FNames=fieldnames(handles.Options);
    FillTo=max(cellNumel(FNames))+2;
    for i=1:numel(FNames)
        Strs{i}=[FNames{i} repmat(' ',[1 FillTo-numel(FNames{i})]) num2str(handles.Options.(FNames{i}))];
    end
    set(handles.Script_Params_listbox,'String',Strs);
    set(handles.Script_Params_listbox,'Value',1);
else
    FNames=setdiff(fieldnames(handles.Options),AdvancedOptions);
    FillTo=max(cellNumel(FNames))+2;
    for i=1:numel(FNames)
        Strs{i}=[FNames{i} repmat(' ',[1 FillTo-numel(FNames{i})]) num2str(handles.Options.(FNames{i}))];
    end
    set(handles.Script_Params_listbox,'String',Strs);
    set(handles.Script_Params_listbox,'Value',1);
end
guidata(hObject, handles);


% --- Executes on selection change in listboxROI.
function listboxROI_Callback(hObject, eventdata, handles)
% hObject    handle to listboxROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxROI contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxROI


% --- Executes during object creation, after setting all properties.
function listboxROI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
CurStr=getKthElement(get(hObject,'String'),get(hObject,'Value'));
if(strcmp(CurStr,'Do nothing'))
    return;
end
if(strcmp(CurStr,'Send report to Gilad'))
    try
        ToWho='giladliberman@gmail.com';
        Subj=['Report for ' handles.ShortSeriesName];
        Msg='Attached.';
        x = inputdlg({'To','Subject','Msg'},...
            'Mail', [1 50;1 50; 5 50],{ToWho,Subj,Msg});
        if(isempty(x))
            return;
        end
        a=gmat2cell(x{3},1);
        for i=1:numel(a);tmp=fliplr(a{i});a{i}=fliplr(tmp(find(tmp~=' ',1):end));end
        gSendMail(x{1},x{2},a,{[handles.destFolder 'Report.pdf']});
    catch
    end
    msgbox('Report sent','Report sent.');
end

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonDestSame.
function pushbuttonDestSame_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDestSame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.destFolderBase=handles.BaseP;

% Save the destination as the default one for next use
base_path_dir = getComputerParams('basepath');
if( exist(base_path_dir,'dir') )
    Dest_Folder_File_Path =[base_path_dir 'Default_Dest_Folder.mat'];
    Default_dest_folder_Val = handles.destFolderBase;
    % Remove old default file
    %if(exist(Dest_Folder_File_Path,'file'))
    %    delete(Dest_Folder_File_Path);
    %end
    save(Dest_Folder_File_Path,'Default_dest_folder_Val');
end
   

set(handles.Dest_Folder_TextBox,'String',handles.destFolderBase);
guidata(hObject, handles);

function DCEMainI=getDCEMainI(handles,GoodS)
Tmp=[handles.ShortInfos(GoodS).ImagesInAcquisition];
Tmp(Tmp==25)=2500; % For PHILIPS now
Tmp(Tmp==55)=2500; % For Verona PHILIPS now
DCEMainI=find(Tmp>50); % For Siemens
if(numel(handles.ShortInfos(1).R1)==1)
    DCEMainI=DCEMainI(1);
end
% If we got more than 1 main, return an error
if(numel(DCEMainI)~=1)
    if(isfield(handles,'ListMode') && handles.ListMode)
        error('Choose DCE Main problem');
    else
        errordlg('Choose DCE Main','DCE main');
    end
    return;
end