
function varargout = MainGUI_Perfusion(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MainGUI_Perfusion_OpeningFcn, ...
    'gui_OutputFcn',  @MainGUI_Perfusion_OutputFcn, ...
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
function MainGUI_Perfusion_OpeningFcn(hObject, eventdata, handles, varargin)
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

% ----------------------------
% Sim_Struct.Use_Model_Selection           = true;
% Sim_Struct.Ignore_Delay_Model_Selection  = false; % Ignore models with delay
% Sim_Struct.Correct_estimation_due_to_delay        = true;       % Try to correct for delay
% Sim_Struct.Min_Time_Delay                         = -0.0;  % Set the minimal possible time delay in seconds for correction
% Sim_Struct.Max_Time_Delay                         = +20.0;  % Set the maximal possible time delay in seconds for correction
% Sim_Struct.Force_RealData_Calc           = true;   % Force the calculation of real data even if calculated bafore
% Sim_Struct.Parallel_Real_Data_Est        = true;   % Parallel the for loop in real data estimation

% Advanced options:

% Sim_Struct.Correct_PVE                   = true;   % Correct Partial Volume Effect for AIF
% Sim_Struct.AIC_Correction                = true; % Use correction for AIC
% Sim_Struct.Data_Weight                   = 0.1;  % Data weight comparing to # of params (Gilad uses 0.1)
% Sim_Struct.poly_deg                      = 4;
% Sim_Struct.knot_interval                 = 5; % knots         = time_vec_minutes(1:knot_interval:end)
% Sim_Struct.LQ_Model_AIF_Delay_Correct             = false;      % Correct AIF by Linear-Quadratic model (Cheong 2003)
% Sim_Struct.Upsampling_resolution_Sec              = 0.05;        % Set the upsampling target
% Sim_Struct.Use_Cyclic_Conv_4_ht_est               = false;       % Use cyclic de-convolution to correct for delay
% Sim_Struct.Simple_AIF_Delay_Correct               = false;       % Correct AIF by max point shift

% ----------------------------

Defaults.Use_Model_Selection                   = 1;
ParamToolTips.Use_Model_Selection              = ['Use_Model_Selection: Run model selection using AIC.' 10 'Default: 1.'];
Defaults.Ignore_Delay_Model_Selection          = 0;
ParamToolTips.Ignore_Delay_Model_Selection     = ['Ignore_Delay_Model_Selection: Run Model Selection without delay as additional parameter.' 10 'Default: 0'];
Defaults.Correct_estimation_due_to_delay       = 1;
ParamToolTips.Correct_estimation_due_to_delay  = ['Correct_estimation_due_to_delay: Use BAT correction in analysis.' 10 'Default: 0'];
Defaults.Min_Time_Delay                        = -0.0;
ParamToolTips.Min_Time_Delay                   = ['Min_Time_Delay: Mix BAT to search.' 10 'Default: -0.0'];
Defaults.Max_Time_Delay                        = +10.0;
ParamToolTips.Max_Time_Delay                   = ['Max_Time_Delay: Max BAT to search.' 10 'Default: +10.0'];
Defaults.Force_RealData_Calc                   = 1;
ParamToolTips.Force_RealData_Calc              = ['Force_RealData_Calc: Force calculation if already calculated.' 10 'Default: 1'];
Defaults.Parallel_Real_Data_Est                = 1;
ParamToolTips.Parallel_Real_Data_Est           = ['Parallel_Real_Data_Est: Run Calculation parallel.' 10 'Default: 1'];
Defaults.Adjusted_Larsson_Model                = 1;
ParamToolTips.Adjusted_Larsson_Model           = ['Use adjusted Larsson Model.' 10 'Default: 1'];

Defaults.Correct_PVE                     = 1;
ParamToolTips.Correct_PVE                = ['Correct_PVE: Correct AIF for PVE.' 10 'Default: 1'];
Defaults.AIC_Correction                  = 1;
ParamToolTips.AIC_Correction             = ['AIC_Correction: Add Correction part to akaike.' 10 'Default: 1'];
Defaults.Data_Weight                     = 0.1;
ParamToolTips.Data_Weight                = ['Data_Weight: Data weight in AICc.' 10 'Default: 0.1'];
Defaults.poly_deg                        = 4;
ParamToolTips.poly_deg                   = ['poly_deg: Degree of polynomial in splines.' 10 'Default: 4'];
Defaults.knot_interval                   = 5;
ParamToolTips.knot_interval              = ['knot_interval: Number of knots in spline calculation' 10 'Default: 5'];
Defaults.LQ_Model_AIF_Delay_Correct      = 0;
ParamToolTips.LQ_Model_AIF_Delay_Correct = ['LQ_Model_AIF_Delay_Correct: Correct BAT using LQ Model.' 10 'Default: 0'];
Defaults.Upsampling_resolution_Sec       = 0.05;
ParamToolTips.Upsampling_resolution_Sec  = ['Upsampling_resolution_Sec: Upsampling resolution for BAT estimation.' 10 'Default: 0.05'];
Defaults.Use_Cyclic_Conv_4_ht_est        = 0;
ParamToolTips.Use_Cyclic_Conv_4_ht_est   = ['Use_Cyclic_Conv_4_ht_est: Use cyclic deconvolution matrix.' 10 'Default: 0'];
Defaults.Simple_AIF_Delay_Correct        = 0;
ParamToolTips.Simple_AIF_Delay_Correct   = ['Simple_AIF_Delay_Correct: Estimate BAT using the simple method.' 10 'Default: 0'];

handles.ParamToolTips = ParamToolTips;

handles.Options = Defaults;
FNames          = fieldnames(handles.Options);
FillTo          = max(cellNumel(FNames))+2;

for i = 1:numel(FNames)
    Strs{i} = [FNames{i} repmat(' ',[1 FillTo-numel(FNames{i})]) num2str(handles.Options.(FNames{i}))];
end

set(handles.Script_Params_Perf_listbox,'String',Strs);
set(handles.Script_Params_Perf_listbox,'Value',1);

tmp=getKthElement(getComputerParams('tpm'),1);
ROIsPath=[tmp(1:end-23) 'Code' filesep 'Utils' filesep 'SPM_precofigures' filesep];

checkboxAdvanced_Perfusion_Callback(hObject, eventdata, handles);

try
    resetGUI(hObject, handles);
catch
end

D=dir([ROIsPath '*.nii']);
ROINames=cellfun(@(x) x(1:end-4),{D.name}','UniformOutput',false);
%set(handles.listboxROI,'String',['Full';ROINames]);
guidata(hObject, handles);

% UIWAIT makes MainGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = MainGUI_Perfusion_OutputFcn(hObject, eventdata, handles)
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

% --- Executes during object creation, after setting all properties.
function Choose_Folder_PushButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Choose_Folder_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in DCE_Analysis_PushButton.
function DCE_Analysis_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to DCE_Analysis_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in Choose_Dest_Folder_PushButton_Perf.
function Choose_Dest_Folder_PushButton_Perf_Callback(hObject, eventdata, handles)
% hObject    handle to Choose_Dest_Folder_PushButton_Perf (see GCBO)
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

% --- Executes on button press in Do_Debug_CheckBox_Perfusion.
function Do_Debug_CheckBox_Perfusion_Callback(hObject, eventdata, handles)
% hObject    handle to Do_Debug_CheckBox_Perfusion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Do_Debug_CheckBox_Perfusion

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

% --- Executes on button press in Run_On_All_Perf_PushButton.
function Run_On_All_Perf_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Run_On_All_Perf_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% ------------------------------------

% Set simulation parameters
Sim_Struct = struct;
Verbosity  = 'None';
Sim_Struct = Simulation_Set_Params(Sim_Struct, Verbosity);

Subject_name          = '';
%Subject_Path          = '\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\DCE-HTR';
Subject_Path          = get(handles.Data_Folder_Text_Box,'String');
PefusionOutput        = get(handles.Dest_Folder_TextBox,'String');

% Override parameters with what user chose
FNames  = {{'Use_Model_Selection'} {'Ignore_Delay_Model_Selection'} {'Correct_estimation_due_to_delay'} ...
           {'Min_Time_Delay'} {'Max_Time_Delay'} {'Force_RealData_Calc'} {'Parallel_Real_Data_Est'} ...
           {'Correct_PVE'} {'AIC_Correction'} {'Data_Weight'} {'Adjusted_Larsson_Model'}...
           {'poly_deg'} {'knot_interval'} {'LQ_Model_AIF_Delay_Correct'} {'Upsampling_resolution_Sec'} ...
           {'Use_Cyclic_Conv_4_ht_est'} {'Simple_AIF_Delay_Correct'}};
       
for i = 1 : length(FNames)
    eval(['Sim_Struct.' cell2mat(FNames{i}) '  = handles.Options.' cell2mat(FNames{i}) ';']);
    %Sim_Struct.(FNames{i})  = handles.Options.(FNames{i});
end

% Run perfusion
Run_DCE_Perfusion( Subject_name, Subject_Path, Sim_Struct, PefusionOutput , Verbosity);

% ------------------------------------



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

function Under_Sampling_EditBox_Perfusion_Callback(hObject, eventdata, handles)
% hObject    handle to Under_Sampling_EditBox_Perfusion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Under_Sampling_EditBox_Perfusion as text
%        str2double(get(hObject,'String')) returns contents of Under_Sampling_EditBox_Perfusion as a double
CurV=get(handles.Script_Params_Perf_listbox,'Value');
Tmp=get(handles.Script_Params_Perf_listbox,'String');
for i=1:numel(Tmp)
    FNames{i}=getKthElement(regexp(Tmp{i},'\W+','split'),1);
end
% FNames=fieldnames(handles.Options);
Str=get(handles.Under_Sampling_EditBox_Perfusion,'String');
Val=str2double(Str);
if(isnan(Val))
    msgbox('Please enter numeric value only.');
    return;
end
i=CurV;
handles.Options.(FNames{i}) = Val;
FillTo=max(cellNumel(FNames))+2;
Strs=get(handles.Script_Params_Perf_listbox,'String');
Strs{i}=[FNames{i} repmat(' ',[1 FillTo-numel(FNames{i})]) num2str(handles.Options.(FNames{i}))];
set(handles.Script_Params_Perf_listbox,'String',Strs);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Under_Sampling_EditBox_Perfusion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Under_Sampling_EditBox_Perfusion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in Script_Params_Perf_listbox.
function Script_Params_Perf_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to Script_Params_Perf_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Script_Params_Perf_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Script_Params_Perf_listbox
CurV = get(handles.Script_Params_Perf_listbox,'Value');
Tmp  = get(handles.Script_Params_Perf_listbox,'String');

% Get the name of the fields
for i = 1:numel(Tmp)
    FNames{i} = getKthElement(regexp(Tmp{i},'\W+','split'),1);
end

% FNames=fieldnames(handles.Options);

% Set the value on handles.options in the relevant field
set(handles.Under_Sampling_EditBox_Perfusion,'String',num2str(handles.Options.(FNames{CurV})));

% Set the value on handles.options in the relevant field
if( numel(CurV) > 0 && isfield( handles.ParamToolTips, FNames{min(CurV)}) )
    set(handles.Script_Params_Perf_listbox,'Tooltip',handles.ParamToolTips.(FNames{min(CurV)}));
end


% --- Executes during object creation, after setting all properties.
function Script_Params_Perf_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Script_Params_Perf_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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

% --- Executes on button press in OpenDest_Perf.
function OpenDest_Perf_Callback(hObject, eventdata, handles)
% hObject    handle to OpenDest_Perf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Explore(handles.destFolderBase);

% --- Executes on button press in checkboxAdvanced_Perfusion.
function checkboxAdvanced_Perfusion_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAdvanced_Perfusion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ----------------------------
% Sim_Struct.Correct_PVE                   = true;   % Correct Partial Volume Effect for AIF
% Sim_Struct.AIC_Correction                = true; % Use correction for AIC
% Sim_Struct.Data_Weight                   = 0.1;  % Data weight comparing to # of params (Gilad uses 0.1)
% Sim_Struct.SNR_single                    = 15; %15 
% Sim_Struct.poly_deg                      = 4;
% Sim_Struct.knot_interval                 = 5; % knots         = time_vec_minutes(1:knot_interval:end)
% Sim_Struct.LQ_Model_AIF_Delay_Correct             = false;      % Correct AIF by Linear-Quadratic model (Cheong 2003)
% Sim_Struct.Upsampling_resolution_Sec              = 0.05;        % Set the upsampling target
% Sim_Struct.Use_Cyclic_Conv_4_ht_est               = false;       % Use cyclic de-convolution to correct for delay
% Sim_Struct.Simple_AIF_Delay_Correct               = false;       % Correct AIF by max point shift
% ----------------------------

% Hint: get(hObject,'Value') returns toggle state of checkboxAdvanced_Perfusion
%AdvancedOptions = {'Correct_PVE','AIC_Correction','Data_Weight','SNR_single','poly_deg','knot_interval','LQ_Model_AIF_Delay_Correct','Upsampling_resolution_Sec','Use_Cyclic_Conv_4_ht_est','Simple_AIF_Delay_Correct'};
AdvancedOptions = {};


if(get(handles.checkboxAdvanced_Perfusion,'Value'))
    FNames = fieldnames(handles.Options);
    FillTo = max(cellNumel(FNames))+2;
    for i=1:numel(FNames)
        Strs{i}=[FNames{i} repmat(' ',[1 FillTo-numel(FNames{i})]) num2str(handles.Options.(FNames{i}))];
    end
    set(handles.Script_Params_Perf_listbox,'String',Strs);
    set(handles.Script_Params_Perf_listbox,'Value',1);
else
    % Get all names in advanced options that do not appear in the regular
    % options
    FNames = setdiff(fieldnames(handles.Options),AdvancedOptions);
    FillTo = max(cellNumel(FNames)) + 2;
    for i=1:numel(FNames)
        Strs{i}=[FNames{i} repmat(' ',[1 FillTo-numel(FNames{i})]) num2str(handles.Options.(FNames{i}))];
    end
    
    set(handles.Script_Params_Perf_listbox,'String',Strs);
    set(handles.Script_Params_Perf_listbox,'Value',1);
    
end
guidata(hObject, handles);



% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
FNames = fieldnames(handles.Options);
FillTo = max(cellNumel(FNames))+2;
for i=1:numel(FNames)
    Strs{i}=[FNames{i} repmat(' ',[1 FillTo-numel(FNames{i})]) num2str(handles.Options.(FNames{i}))];
end
set(handles.Script_Params_Perf_listbox,'String',Strs);
set(handles.Script_Params_Perf_listbox,'Value',1);


