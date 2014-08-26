function varargout = bolus_properties_GUI(varargin)
% BOLUS_PROPERTIES_GUI MATLAB code for bolus_properties_GUI.fig
%      BOLUS_PROPERTIES_GUI, by itself, creates a new BOLUS_PROPERTIES_GUI or raises the existing
%      singleton*.
%
%      H = BOLUS_PROPERTIES_GUI returns the handle to a new BOLUS_PROPERTIES_GUI or the handle to
%      the existing singleton*.
%
%      BOLUS_PROPERTIES_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BOLUS_PROPERTIES_GUI.M with the given input arguments.
%
%      BOLUS_PROPERTIES_GUI('Property','Value',...) creates a new BOLUS_PROPERTIES_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bolus_properties_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bolus_properties_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bolus_properties_GUI

% Last Modified by GUIDE v2.5 01-May-2014 11:54:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bolus_properties_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @bolus_properties_GUI_OutputFcn, ...
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


% --- Executes just before bolus_properties_GUI is made visible.
function bolus_properties_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bolus_properties_GUI (see VARARGIN)

% Choose default command line output for bolus_properties_GUI
handles.output = hObject;

% set(handles.figure1,'tag','bolus_properties_GUI');
axes(handles.DSC_curve_fig);
curr_folder=pwd;
imshow([curr_folder,'\DSC_curve.JPG']);
first_bl_sample=getappdata(0,'first_bl_sample');
set(handles.first_bl_sample_value,'String',num2str(first_bl_sample));
last_bl_sample=getappdata(0,'last_bl_sample');
set(handles.last_bl_sample_value,'String',num2str(last_bl_sample));
last_sample=getappdata(0,'last_sample');
set(handles.last_sample_value,'String',num2str(last_sample));



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bolus_properties_GUI wait for user response (see UIRESUME)
% uiwait(handles.bolus_properties_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = bolus_properties_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function first_bl_sample_value_Callback(hObject, eventdata, handles)
% hObject    handle to first_bl_sample_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of first_bl_sample_value as text
%        str2double(get(hObject,'String')) returns contents of first_bl_sample_value as a double

% --- Executes during object creation, after setting all properties.
function first_bl_sample_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to first_bl_sample_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function last_bl_sample_value_Callback(hObject, eventdata, handles)
% hObject    handle to last_bl_sample_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of last_bl_sample_value as text
%        str2double(get(hObject,'String')) returns contents of last_bl_sample_value as a double


% --- Executes during object creation, after setting all properties.
function last_bl_sample_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to last_bl_sample_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function last_sample_value_Callback(hObject, eventdata, handles)
% hObject    handle to last_sample_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of last_sample_value as text
%        str2double(get(hObject,'String')) returns contents of last_sample_value as a double


% --- Executes during object creation, after setting all properties.
function last_sample_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to last_sample_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ok_button.
function ok_button_Callback(hObject, eventdata, handles)
% hObject    handle to ok_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.user_pushed_button=1;

first_bl_sample=str2double(get(handles.first_bl_sample_value,'String'));
last_bl_sample=str2double(get(handles.last_bl_sample_value,'String'));
last_sample=str2double(get(handles.last_sample_value,'String'));

hMainGui=getappdata(0,'hMainGui');
setappdata(hMainGui,'first_bl_sample',first_bl_sample);
setappdata(hMainGui,'last_bl_sample',last_bl_sample);
setappdata(hMainGui,'last_sample',last_sample);

guidata(hObject, handles);


close(handles.bolus_properties_GUI) ;

% --- Executes during object creation, after setting all properties.
function DSC_curve_fig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DSC_curve_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate DSC_curve_fig


% --- Executes during object creation, after setting all properties.
function DSC_curve_fig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DSC_curve_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate DSC_curve_fig
