
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

% If using FMRI-T9, default output folder will be the desktop
if strcmp('FMRI-T9',getenv('ComputerName'))
        handles.destFolderBase= 'C:\Users\guyn\Desktop\DCE_OUT\';
        % Set the Temp directory under the destination directory
        setComputerParamM('temppath',[handles.destFolderBase 'Temp']);
end

set(handles.Dest_Folder_TextBox,'String',handles.destFolderBase);
set(handles.list_of_paths_box,'Min',0);
set(handles.list_of_paths_box,'Max',3);

Defaults.SubSampling=1;
Defaults.nVolsToRemoveFromEnd=0;
Defaults.SubSecondResolution=2;
Defaults.MinFirstBolusStd=2;
Defaults.EM_Num_Of_Iterations=5;
Defaults.FMS_TolFun=1e-11;
Defaults.FMS_MaxFunEvals=10000;
Defaults.FMS_MaxIter=10000;
Defaults.MaxTDif_ForAIFSearch=3;
Defaults.MaxTDif_ForWholeVOI=6;
Defaults.Rep_MaxAroundBolus=10;
Defaults.Rep_RatioToEnd=10;
Defaults.Rep_nPerSet=1;
Defaults.MakeNoBATManualArtAnalysis=0;
Defaults.MakeBATManualArtAnalysis=0;
Defaults.MakeBATAutoArtAnalysis=1;

handles.Options=Defaults;
FNames=fieldnames(handles.Options);
FillTo=max(cellNumel(FNames))+2;
for i=1:numel(FNames)
    Strs{i}=[FNames{i} repmat(' ',[1 FillTo-numel(FNames{i})]) num2str(handles.Options.(FNames{i}))];
end
set(handles.Script_Params_listbox,'String',Strs);
set(handles.Script_Params_listbox,'Value',1);

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
if(isfield(handles,'BaseP'))
    handles.BaseP=uigetdir(handles.BaseP);
else

    % Default directory when Guy uses the GUI
    % At Home
    if strcmp('GUYOV-PC',getenv('ComputerName'))
        handles.BaseP=uigetdir('H:\Guy\Dropbox\University\Msc\Thesis\Matlab\Data');

        % At University
    elseif strcmp('FMRI-T9',getenv('ComputerName'))
        %handles.BaseP=uigetdir('\\FMRI-GUY2\Dropbox\University\Msc\Thesis\SourceForge\Development\Data');
        handles.BaseP=uigetdir('Z:\WM-lesions\BAVLI_YEHOSHUA\Study20121213_084054\');

        % Default
    else
        % At Ichilov - Unix
        handles.BaseP=uigetdir('/data/');
    end

end

handles.BaseP=[handles.BaseP filesep];
set(handles.Data_Folder_Text_Box,'String',handles.BaseP);

MatFN=[handles.BaseP 'Data.mat'];
if(exist(MatFN,'file'))
    handles=Read_Protocol_Action_Callback(hObject, eventdata, handles);
end

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


% --- Executes on button press in Process4D_PushButton.
function handles=Process4D_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to Process4D_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function handles = Process4D(hObject, eventdata, handles)

% Get indices of the data in the listbox
GoodS=get(handles.DCE_Data_ListBox,'Value');

% The following is used to get the DCE over time data
DCEMainI=find([handles.ShortInfos(GoodS).ImagesInAcquisition]>100);
% If we got more than 1 main, return an error
if(numel(DCEMainI)~=1)
    if(isfield(handles,'ListMode') && handles.ListMode)
        error('Choose DCE Main problem');
    else
        errordlg('Choose DCE Main','DCE main');
    end
    return;
end

f = helpdlg('Preprocess...');

% Create relevant masks for the 4D data and find bolus time.
DCET1_Prepare4Df(handles.ShortInfos(GoodS(DCEMainI)),handles.destFolder,get(handles.forceP4d,'Value'),handles.Options);

% If help dialog is still open, close it
if ( exist('f','var') )
    close(f);
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
DCEMainI=find([handles.ShortInfos(GoodS).ImagesInAcquisition]>100);
% Remove the Main from the rest
GoodS2=GoodS(setdiff(1:numel(GoodS),DCEMainI));
% Get the Info of the relaxometry data
RelaxInfox=handles.ShortInfos(GoodS2(ismember({handles.ShortInfos(GoodS2).Class},{'DCE','T1MAP'})));
% Check of N3 was chosen
DoN3=get(handles.N3_B1_CheckBox,'Value');
% ASK GILAD - what is the purpose of DoGlobal and DoDeviations?
% ANSWER - 
% DoDeviations � Correction for flip angles (exp. � machine says its 20 degrees while the real one is 20.2). 
%                                     Implemented under CalcRelaxForVolNFA.m. Currently not used.
%
% DoGlobal � N3 fixes B1 non-uniformity by a factor for each voxel. 
%                          Gilad picks a point in white matter to be the reference point (implemented in DCET1_RelaxForSubjectf.m).
DoGlobal=get(handles.Global_FA_Factor_CheckBox,'Value');
DoDeviations=get(handles.FA_Deviations_CheckBox,'Value');
% See if we chose to force calc T1 (if not, calculation is skipped if was previously done)
CalcForce=get(handles.forcet1,'Value');

% Get from GUI to what to co-register T1 maps 
CrgToS=get(handles.t1coregto,'String');
CrgTo=CrgToS{get(handles.t1coregto,'Value')};

% Starting T1 calculation
f = helpdlg('T1 calculation...');

% Calculate T1 maps
DCET1_RelaxForSubjectf(RelaxInfox,handles.destFolder,DoN3,DoGlobal,DoDeviations,CalcForce,CrgTo,handles.Options);

% If help dialog is still open, close it
if ( exist('f','var') )
    close(f);
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
DoN3=get(handles.N3_B1_CheckBox,'Value');
DoGlobal=get(handles.Global_FA_Factor_CheckBox,'Value');
DoDeviations=get(handles.FA_Deviations_CheckBox,'Value');
CalcForce=get(handles.forcectc,'Value');

GoodS=get(handles.DCE_Data_ListBox,'Value');
DCEMainI=find([handles.ShortInfos(GoodS).ImagesInAcquisition]>100);

% UnderSampling=str2num(get(handles.Under_Sampling_EditBox,'String'));
UnderSampling=handles.Options.SubSampling;

f = helpdlg('CTCs calculation...');
DCET1_CTCf(handles.ShortInfos(GoodS(DCEMainI)),handles.destFolder,DoN3,DoGlobal,DoDeviations,CalcForce,handles.Options);
% If help dialog is still open, close it
if ( exist('f','var') )
    close(f);
end

% --- Executes on button press in DCE_Analysis_PushButton.
function DCE_Analysis_PushButton_Callback(hObject, eventdata, handles)
% hObject    handle to DCE_Analysis_PushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
DCEMainI=find([handles.ShortInfos(GoodS).ImagesInAcquisition]>100);
CurMainDCEInfo=handles.ShortInfos(GoodS(DCEMainI));
f = helpdlg('AIF extraction...');
disp('DCET1_PKf..');
DCET1_PKf;
disp('AIFFinderTry..');
AIFFinderTry2;
% If help dialog is still open, close it
if ( exist('f','var') )
    close(f);
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
DCEMainI=find([handles.ShortInfos(GoodS).ImagesInAcquisition]>100);
CurMainDCEInfo=handles.ShortInfos(GoodS(DCEMainI));

f = helpdlg('Whole ROI computation...');
disp('DCET1_PKf..');
% Get clustering, bolus time and noise...
DCET1_PKf;
disp('WholeROICompute..');
%
DCET1_WholeVolCompute;
% If help dialog is still open, close it
if ( exist('f','var') )
    close(f);
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

MatFN=[CurBaseP 'Data.mat'];

% Pop help dialog
f = helpdlg('Loading...');

% ASK GILAD - when I use "Force Read Protocol" I get a failure in that condition
%
if(exist(MatFN,'file') && ~get(handles.Force_Read_Prot_CheckBox,'Value')) %% If exists and not "force"
    load(MatFN);
    handles.ShortInfos=ShortInfos;
else % Either "force" or .mat file does not exist

    if(~exist('ShortInfos','var'))
        AddNewScans(handles.BaseP,1);
        List_From_DataBase_PushButton_Callback(hObject, eventdata, handles);
    end

    ShortInfos=getSpecificInfos('xPath',CurBaseP);

    if(isempty(ShortInfos))
        ShortInfos = ReadNewScans(CurBaseP);
    end
    handles.ShortInfos=ShortInfos;
    save(MatFN,'ShortInfos');
end

% If help dialog is still open, close it
if ( exist('f','var') )
    close(f);
end

set(handles.Data_Folder_Text_Box,'String',[handles.ShortInfos(end).Name ' ' handles.ShortInfos(end).SeriesDate]);
Descriptions=strcat({handles.ShortInfos.SeriesDescription},'_');
Descriptions=strcat(Descriptions,{handles.ShortInfos.Class});
Descriptions=strrep(Descriptions,'_',' ');

DCEMainI=find([handles.ShortInfos.ImagesInAcquisition]>100 & strcmp({handles.ShortInfos.Class},'DCE'));
for i=DCEMainI
    Descriptions{i}=[Descriptions{i} ' ' num2str(handles.ShortInfos(i).FlipAngle) ' ' num2str(handles.ShortInfos(i).ImagesInAcquisition)];
end

set(handles.DCE_Data_ListBox,'String',Descriptions);
set(handles.DCE_Data_ListBox,'Max',3);
set(handles.DCE_Data_ListBox,'Min',0);

% Shorten the series name
handles.ShortSeriesName=[ToShortName(handles.ShortInfos(end).Name) '_' handles.ShortInfos(end).SeriesDate];
% Remove any spaces
handles.ShortSeriesName = regexprep(handles.ShortSeriesName,'\W','')
% Change the destination directory to the short name
handles.destFolder=[handles.destFolderBase handles.ShortSeriesName filesep];

set(handles.Dest_Folder_TextBox,'String',handles.destFolder);
[s, mess, messid] =mkdir(handles.destFolder);

handles.GoodSers=find(ismember({handles.ShortInfos.Class},{'T1MAP','DCE'})); % 'FIESTA'

% GUYN - DATA CHANGE - commented to select all T1 data and not just up to Main data

% % Get the index of the DCE Main data
% DCEMainI=find([handles.ShortInfos(handles.GoodSers).ImagesInAcquisition]>100);
% % If we have DCE Main acquisition, mark all data up to it.
% if(~isempty(DCEMainI))
%     handles.GoodSers=[handles.GoodSers(1:(DCEMainI(1)-1)) handles.GoodSers(DCEMainI(end))];
% end

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
Strs=get(handles.Stages_List_Box,'String');
V=get(handles.Stages_List_Box,'Value');
ToDo=Strs(V);
if(~exist(handles.destFolderBase,'dir'))
    mkdir(handles.destFolderBase);
end

if(ismember('Process 4D',ToDo))
    handles=Process4D(hObject, eventdata, handles);
end
if(ismember('Process T1',ToDo))
    handles=ProcessT1(hObject, eventdata, handles);
end
if(ismember('Create CTCs',ToDo))
    handles=CreateCTCs(hObject, eventdata, handles);
end
if(ismember('DCE Analysis',ToDo))
    AIFFind(hObject, eventdata, handles);
end
if(ismember('Whole ROI computation',ToDo))
    WholeROICompute(hObject, eventdata, handles);
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
DoDebug=get(handles.Do_Debug_CheckBox,'Value');
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
    DDC3{i}=DDC2{i}(1:S(end-1));
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
DCEMainI=find([handles.ShortInfos(GoodS).ImagesInAcquisition]>100);
CurMainDCEInfo=handles.ShortInfos(GoodS(DCEMainI));
disp('DCET1_PKf..');
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

AA=getSpecificInfos('Class','DCE');
BB=AA([AA.ImagesInAcquisition]>100);
for i=1:numel(BB)
    Tmp=find(BB(i).Path==filesep);
    Before{i}=BB(i).Path(1:Tmp(end));
end
CurPs=unique(Before)';
set(handles.list_of_paths_box,'String',CurPs);
handles.CurPs=CurPs;
% If help dialog is still open, close it
if ( exist('f','var') )
    close(f);
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

    DCEMainI=find([ShortInfos.ImagesInAcquisition]>100 & strcmp({ShortInfos.Class},'DCE'));
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
if ( exist('f','var') )
    close(f);
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
save(handles.GUISaveFN,'DestFolderBase');
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
FNames=fieldnames(handles.Options);
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
FNames=fieldnames(handles.Options);
set(handles.Under_Sampling_EditBox,'String',num2str(handles.Options.(FNames{CurV})));


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


