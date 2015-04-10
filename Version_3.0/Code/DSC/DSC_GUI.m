function varargout = DSC_GUI(varargin)
% DSC_GUI MATLAB code for DSC_GUI.fig
%      DSC_GUI, by itself, creates a new DSC_GUI or raises the existing
%      singleton*.
%
%      H = DSC_GUI returns the handle to a new DSC_GUI or the handle to
%      the existing singleton*.
%
%      DSC_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSC_GUI.M with the given input arguments.
%
%      DSC_GUI('Property','Value',...) creates a new DSC_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DSC_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DSC_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DSC_GUI

% Last Modified by GUIDE v2.5 14-Aug-2014 16:56:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DSC_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DSC_GUI_OutputFcn, ...
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

% --- Executes just before DSC_GUI is made visible.
function DSC_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DSC_GUI (see VARARGIN)

% Choose default command line output for DSC_GUI
handles.output = hObject;

axes(handles.DSC_curve_fig);
imshow('D:\users\chenh\DSC_project\development\DSC_curve.JPG')

% axes(handles.data_window);
time_sample=get(handles.time_sample_slider,'Value');
slice=get(handles.slice_slider,'Value');
clc;

% Update the status of parameters checkboxes for SVD:
if get(handles.checkbox_sSVD,'Value')
    set(handles.checkbox_sSVD_th,'Enable','on');
    set(handles.text_sSVD_th,'Enable','on');
end
if get(handles.checkbox_cSVD,'Value')
    set(handles.checkbox_cSVD_th,'Enable','on');
    set(handles.text_cSVD_th,'Enable','on');
end
if get(handles.checkbox_oSVD,'Value')
    set(handles.checkbox_oSVD_OI,'Enable','on');
    set(handles.text_oSVD_OI,'Enable','on');
end

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes DSC_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DSC_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_CoRegistration.
function checkbox_CoRegistration_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_CoRegistration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_CoRegistration


% --- Executes on button press in checkbox_LowAtZero.
function checkbox_LowAtZero_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_LowAtZero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.filter_LowAtZero_enable=get(hObject,'Value');
guidata(hObject,handles);


% --- Executes on button press in checkbox_LowAverageBaseline.
function checkbox_LowAverageBaseline_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_LowAverageBaseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.filter_LowAverageBaseline_enable=get(hObject,'Value');
guidata(hObject,handles);


% --- Executes on button press in checkbox_NoBolusPeak.
function checkbox_NoBolusPeak_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_NoBolusPeak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.filter_NoBolusPeak_enable=get(hObject,'Value');
guidata(hObject,handles);


% --- Executes on button press in checkbox_PeakOutOfBolusRange.
function checkbox_PeakOutOfBolusRange_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PeakOutOfBolusRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.filter_PeakOutOfBolusRange_enable=get(hObject,'Value');
guidata(hObject,handles);


% --- Executes on button press in checkbox_BigFluctuations.
function checkbox_BigFluctuations_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_BigFluctuations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.filter_BigFluctuations_enable=get(hObject,'Value');
guidata(hObject,handles);


% --- Executes on button press in pre_process.
function pre_process_Callback(hObject, eventdata, handles)
% hObject    handle to pre_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guidata(hObject);
%Display message:
set(handles.processing_status_display,'String','Processing . . .');
% Filters for pre-process
% checkbox_LowAtZero_Callback(hObject, eventdata, handles);
% checkbox_LowAverageBaseline_Callback(hObject, eventdata, handles);
% checkbox_NoBolusPeak_Callback(hObject, eventdata, handles);
% checkbox_PeakOutOfBolusRange_Callback(hObject, eventdata, handles);
% checkbox_BigFluctuations_Callback(hObject, eventdata, handles);
handles=guidata(hObject);
enable2=get(handles.checkbox_LowAtZero,'Value');%    handles.filter_LowAtZero_enable;
enable3=get(handles.checkbox_LowAverageBaseline,'Value');
enable4=get(handles.checkbox_NoBolusPeak,'Value');
enable5=get(handles.checkbox_PeakOutOfBolusRange,'Value');
enable6=get(handles.checkbox_BigFluctuations,'Value');
enable7=get(handles.checkbox_ZeroValues,'Value');
filter_is_on=[enable2 enable3 enable4 enable5 enable6 enable7];

% Apply the filters
data4D_init=handles.data4D_init;
Header=handles.header;
[Mask,data4D_filtered,Filter_per_type,Filter_total] = data4D_mask_filtering(data4D_init,filter_is_on);
num_voxels_total=numel(Mask);
data_voxels_row_inds=Filter_total.data_voxels_row_col_slice(:,1);
data_voxels_col_inds=Filter_total.data_voxels_row_col_slice(:,2);
data_voxels_slice_inds=Filter_total.data_voxels_row_col_slice(:,3);
handles.data4D_filtered=data4D_filtered;
handles.Mask=Mask;
handles.data_voxels_row_col_slice=Filter_total.data_voxels_row_col_slice;
handles.data_voxels_indices=Filter_total.data_voxels_indices;
guidata(hObject,handles);

%Show several DSC data time curves on a separate figure:
num_DSC_curves_to_show=30;
DSC_curve_indices_to_show=ceil(length(Filter_total.data_voxels_indices)*rand([1 num_DSC_curves_to_show]));
DSC_curves_fig=figure;
for ii=1:num_DSC_curves_to_show
   plot(squeeze(log(data4D_init(data_voxels_row_inds(DSC_curve_indices_to_show(ii)),data_voxels_col_inds(DSC_curve_indices_to_show(ii)),data_voxels_slice_inds(DSC_curve_indices_to_show(ii)),:))));
   hold all;
end


% convert the filtered 4D data to concentration c(t):   %%%%%%%%%%%%
name_str=Header.hdr.hk.db_name;
TE=str2double(name_str(strfind(name_str,'TE')+3:strfind(name_str,'TE')+4));
if TE<0 || TE>1000
    error('TE has an invalid value: ',TE);
end
update_bolus_properties_Callback(hObject, eventdata, handles);
handles=guidata(hObject);
baseline_edges=[handles.first_bl_sample handles.last_bl_sample];
[concentration_4D Mask] = Intens2concentration_4D(data4D_filtered,baseline_edges,TE,Mask);

% find global min/max concentrations for later use: (only from the voxels
% that passed the mask)
majority_perc=90; % choose % of max voxels to be below a certain level to be shown in the GUI
concentration_4D_vec=concentration_4D(:);
concent_4D_vec_sorted=sort(concentration_4D_vec,'descend');
concent_4D_vec_sorted_non_zero=concent_4D_vec_sorted(find(concent_4D_vec_sorted>0));
ct_small_than_inf_inds=find(concent_4D_vec_sorted_non_zero<Inf);
c_t_max_tot=concent_4D_vec_sorted_non_zero(1);%concent_4D_vec_sorted(ct_small_than_inf_inds(1));
c_t_max=concent_4D_vec_sorted_non_zero(ct_small_than_inf_inds(ceil((1-majority_perc/100)*length(Filter_total.data_voxels_indices))));
ct_large_than_minus_inf_inds=find(concent_4D_vec_sorted>-Inf);
c_t_min=concent_4D_vec_sorted(ct_large_than_minus_inf_inds(end));
handles.c_t_global_min=c_t_min;
handles.c_t_global_max=c_t_max_tot;
handles.c_t_majority_max=c_t_max;
handles.curves4aif_num=0;
handles.inds4aif_mat=zeros(1,3);
handles.AIF=zeros(1,handles.N_time_points_final);

% plot initial data in the data window section
data4D_rgb=repmat(data4D_init,[1 1 1 1 3]);
axes(handles.data_window);
time_sample=get(handles.time_sample_slider,'Value');
slice=get(handles.slice_slider,'Value');
% handles.slice_image_obj=image(mat2gray(data4D_init(:,:,slice,time_sample)));
image_to_show=squeeze(data4D_rgb(:,:,slice,time_sample,:));
imshow(mat2gray(image_to_show));
% set(handles.slice_image_obj,'ButtonDownFcn', {@slice_image_ButtonDownFcn});
set(gcf,'WindowButtonMotionFcn', {@mouseMoveOverDataWindow,handles.data_window});
set(gcf,'WindowButtonDownFcn', {@mousePressDataWindow,handles.data_window});

mouse_pos=get(handles.data_window,'currentPoint');
col=ceil(mouse_pos(1,1)); %col
row=ceil(mouse_pos(1,2)); %row
if 0<col && col<=size(data4D_init,1) && 0<row && row<=size(data4D_init,2)
    set(handles.row_text,'String',num2str(row));
    set(handles.col_text,'String',num2str(col));
end
% handles=guidata(hObject);
% pos=handles.mouse_pos;
% data_windows_position=get(handles.data_window,'position');
% handles.pixelx_to_show_ct=ceil(size(data4D_init,1)*(pos(1)-data_windows_position(1))/(data_windows_position(3)));
% handles.pixely_to_show_ct=ceil(size(data4D_init,2)*(pos(2)-data_windows_position(2))/(data_windows_position(4)));
% set(handles.row_text,'String',num2str(handles.pixelx_to_show_ct));
% set(handles.col_text,'String',num2str(handles.pixely_to_show_ct));
% disp(['You clicked X:',num2str(pos(1)),', Y:',num2str(pos(2))]);

% save data to handles:
handles.concentration_4D=concentration_4D;
handles.time_sample=time_sample;
handles.data4D_rgb=data4D_rgb;

guidata(hObject,handles);

% display "done" message:
set(handles.processing_status_display,'String',['Done. Num of voxels passed mask: ',num2str(uint8(100*length(Filter_total.data_voxels_indices)/num_voxels_total)) ,'% .']);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in nii_button.
function nii_button_Callback(hObject, eventdata, handles)
% hObject    handle to nii_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val=get(hObject,'value');
if val==get(hObject,'Max')
    handles.data_file_type='nii';
else
    handles.data_file_type='dicom';
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function OpenDataChooseType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OpenDataChooseType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in opendata.
function opendata_Callback(hObject, eventdata, handles)
% hObject    handle to opendata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nii_button_Callback(hObject, eventdata, handles);
handles=guidata(hObject);
file_type=handles.data_file_type;
if strcmp(file_type,'dicom')
    file_type_ext='dcm';
else 
    file_type_ext=file_type;
end
file_type_str=strcat('*.',file_type_ext);
[FileName,PathName] = uigetfile({file_type_str});
data_file_path=strcat(PathName,FileName);
set(handles.file_name_to_display,'String',strcat('Opened file: ',data_file_path));
if strcmp(handles.data_file_type,'nii')
    [Out Header A]=loadniidata(data_file_path);
else
    % read Dicom is not support yet
    Display('Reading other file than nii is not supported yet');
end
% arrange the image data in the proper orientation:
voxels_data_4D=flipdim(permute(Out,[2 1 3 4]),1);

handles.data4D_init=voxels_data_4D;
handles.header=Header;
guidata(hObject, handles);
% Set the values of the data window sliders according to the data:
N_slices=size(voxels_data_4D,3);
N_time_points_init=size(voxels_data_4D,4);
set(handles.slice_slider,'Enable','on');
set(handles.slice_slider,'Max',N_slices);
set(handles.slice_slider,'Min',1);
set(handles.slice_slider,'SliderStep',[1 1]/(N_slices-1));
set(handles.slice_slider,'Value',1);
set(handles.time_sample_slider,'Enable','on');
set(handles.time_sample_slider,'Max',N_time_points_init);
set(handles.time_sample_slider,'Min',1);
set(handles.time_sample_slider,'SliderStep',[1 1]/(N_time_points_init-1));
set(handles.time_sample_slider,'Value',1);

guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function file_name_to_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_name_to_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



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


% --- Executes during object creation, after setting all properties.
function DSC_curve_fig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DSC_curve_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate DSC_curve_fig


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Manual_AIF_button.
function Manual_AIF_button_Callback(hObject, eventdata, handles)
% hObject    handle to Manual_AIF_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in update_bolus_properties.
function update_bolus_properties_Callback(hObject, eventdata, handles)
% hObject    handle to update_bolus_properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

% Calc C(t) curves with the updated bolus properties:
first_bl_sample=str2double(get(handles.first_bl_sample_value,'String'));
last_bl_sample=str2double(get(handles.last_bl_sample_value,'String'));
last_sample=str2double(get(handles.last_sample_value,'String'));
baseline_edges=[first_bl_sample last_bl_sample];
Header=handles.header;
Mask=handles.Mask;
name_str=Header.hdr.hk.db_name;
TE=str2double(name_str(strfind(name_str,'TE')+3:strfind(name_str,'TE')+4));
if TE<0 || TE>1000
    error('TE has an invalid value: ',TE);
end
data4D_filtered=handles.data4D_filtered;
[concentration_4D] = Intens2concentration_4D(data4D_filtered,baseline_edges,TE,Mask);

handles.concentration_4D=concentration_4D;
handles.first_bl_sample=first_bl_sample;
handles.last_bl_sample=last_bl_sample;
handles.last_sample=last_sample;
handles.N_time_points_final=last_sample-first_bl_sample+1;

%Update time samples slider of data window:
max_slider_value=last_sample;
min_slider_value=first_bl_sample;
set(handles.time_sample_slider,'Value',min_slider_value);
set(handles.time_sample_slider,'Max',max_slider_value);
set(handles.time_sample_slider,'Min',min_slider_value);
set(handles.time_sample_slider,'SliderStep',[1 1]/(max_slider_value-min_slider_value));

guidata(hObject,handles);


% --- Executes on button press in open_data_window.
function open_data_window_Callback(hObject, eventdata, handles)
% hObject    handle to open_data_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on slider movement.
function slice_slider_Callback(hObject, eventdata, handles)
% hObject    handle to slice_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles=guidata(hObject);
slice=round(get(hObject,'Value'));
set(handles.slice_text,'String',num2str(slice));
axes(handles.data_window);
time_sample=round(get(handles.time_sample_slider,'Value'));
data4D_rgb=handles.data4D_rgb;
image_to_show=squeeze(data4D_rgb(:,:,slice,time_sample,:));
imshow(mat2gray(image_to_show));


% --- Executes during object creation, after setting all properties.
function slice_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slice_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function time_sample_slider_Callback(hObject, eventdata, handles)
% hObject    handle to time_sample_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles=guidata(hObject);
time_sample=round(get(hObject,'Value'));
set(handles.timepoint_text,'String',num2str(time_sample));
axes(handles.data_window);
slice=round(get(handles.slice_slider,'Value'));
data4D_rgb=handles.data4D_rgb;
image_to_show=squeeze(data4D_rgb(:,:,slice,time_sample,:));
imshow(mat2gray(image_to_show));

handles.time_sample=time_sample;
guidata(hObject);


% --- Executes during object creation, after setting all properties.
function time_sample_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_sample_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% N_time_points_to_show=str2double(get(handles.last_sample_value,'String'))-str2double(get(handles.first_bl_sample_value,'String'))+1;

% --- Executes during object creation, after setting all properties.
function slice_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slice_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function timepoint_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timepoint_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function data_window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate data_window
% set(handles.data_window,'WindowButtonMotionFcn', @mouseMoveOverDataWindow);


% --- Executes during object creation, after setting all properties.
function processing_status_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to processing_status_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function row_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to row_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function col_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to col_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function mouseMoveOverDataWindow(hObject,eventdata,handles)

handles=guidata(hObject);
pos = get(handles.data_window,'CurrentPoint');
% disp(['Your mouse is on X:',num2str(pos(1,1)),', Y:',num2str(pos(1,2))]);
% handles.mouse_pos=pos;
% hObject.mouse_pos=pos;
col=ceil(pos(1,1)); %col
row=ceil(pos(1,2)); %row
if 0<col && col<=size(handles.data4D_init,1) && 0<row && row<=size(handles.data4D_init,2)
    set(handles.row_text,'String',num2str(row));
    set(handles.col_text,'String',num2str(col));
    axes(handles.plot_AIFs);
    slice=round(get(handles.slice_slider,'Value'));
    concentration_4D=handles.concentration_4D;
    Mask=handles.Mask;
    AIF=handles.AIF;
    c_t_curve_to_show=squeeze(concentration_4D(row,col,slice,:));
    % choose the max of y axis, for a comfortable view (there are 2
    % different options. One is the global maximum, the other is the max of
    % ~90% from the curves.
    local_max=max(c_t_curve_to_show);
    if local_max>handles.c_t_majority_max
        y_axis_max=handles.c_t_global_max;
    else
        y_axis_max=handles.c_t_majority_max;
    end
    plot(AIF,'r','LineWidth',1.5);
    hold on;
    plot(c_t_curve_to_show);
    hold off;
    axis([handles.first_bl_sample handles.last_sample handles.c_t_global_min y_axis_max]);
%     axis([0 length(concentration_4D(row,col,slice,:)) c_t_min-5 c_t_max+5]);
end

function mousePressDataWindow(hObject,eventdata,handles)

handles=guidata(hObject);
data4D=handles.data4D_init;
time_sample=handles.time_sample;
% get coordinates:
pos = get(handles.data_window,'CurrentPoint');
curves4aif_num=handles.curves4aif_num;
inds4aif_mat=handles.inds4aif_mat;
col=ceil(pos(1,1)); %col
row=ceil(pos(1,2)); %row
% if mouse is inside the image:
%     - choose the curve of the voxel chosen.
if 0<col && col<=size(handles.data4D_init,1) && 0<row && row<=size(handles.data4D_init,2)

    axes(handles.plot_AIFs);
    slice=round(get(handles.slice_slider,'Value'));
    concentration_4D=handles.concentration_4D;
    Mask=handles.Mask;
    curves4aif_num=curves4aif_num+1;  % update the number of curves participating
    inds4aif_mat(curves4aif_num,:)=[row col slice];
    % compute and plot AIF with the new concentraion curve:
    [AIF,AIF_std]=calc_AIF_from_curves(concentration_4D,inds4aif_mat,'average');
    % choose the max of y axis, for a comfortable view (there are 2
    % different options. One is the global maximum, the other is the max of
    % ~90% from the curves.
    if max(AIF)>handles.c_t_majority_max
        y_axis_max=handles.c_t_global_max;
    else
        y_axis_max=handles.c_t_majority_max;
    end
    plot(AIF,'r','Linewidth',1.5);
%     hold on;
    axis([handles.first_bl_sample handles.last_sample handles.c_t_global_min y_axis_max]);    
    
   % mark the chosen voxel in the brain image (data_window)
   
    axes(handles.data_window);
    data4D_rgb=handles.data4D_rgb;
    %mark voxel in yellow:
    image_max=max(max(max(data4D_rgb(:,:,slice,time_sample,:))));
    yellow_value=image_max*reshape([1 1 0],1,1,1,1,3);
    data4D_rgb(row,col,slice,time_sample,:)=yellow_value;
    image_to_show=squeeze(data4D_rgb(:,:,slice,time_sample,:));
    imshow(mat2gray(image_to_show));

    % update curves listbox:
    String_chen={'Row 5';'Row 5'};
    curves4aif_listbox_string=cell(curves4aif_num,1);
    for ii=1:size(inds4aif_mat,1);
        curves4aif_listbox_string{ii}=['Row ',num2str(inds4aif_mat(ii,1)),'    Col ',num2str(inds4aif_mat(ii,2)),'    Slice ',num2str(inds4aif_mat(ii,3))];
        
    end
    set(handles.curves4aif_listbox,'String',curves4aif_listbox_string);
    
    % update values in handles:
    handles.curves4aif_num=curves4aif_num;
    handles.inds4aif_mat=inds4aif_mat;
    handles.AIF=AIF;
    handles.data4D_rgb=data4D_rgb;
end
guidata(hObject,handles);


% --- Executes on mouse press over axes background.
function data_window_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to data_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
pos = get(hObject,'CurrentPoint');
hObject.mouse_pos=pos;
guidata(hObject,handles);


function slice_image_ButtonDownFcn(hObject, eventdata, handles)
handles=guidata(hObject);
pos = get(hObject,'CurrentPoint');
hObject.mouse_pos=pos;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function plot_AIFs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_AIFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate plot_AIFs


% --- Executes on selection change in curves4aif_listbox.
function curves4aif_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to curves4aif_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns curves4aif_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from curves4aif_listbox


% --- Executes during object creation, after setting all properties.
function curves4aif_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to curves4aif_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in produce_maps_button.
function produce_maps_button_Callback(hObject, eventdata, handles)
% hObject    handle to produce_maps_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
AIF=handles.AIF;
concentration4D=handles.concentration_4D;
Mask=handles.Mask;
deconv_methods.sSVD.en=get(handles.checkbox_sSVD,'Value');
deconv_methods.sSVD.th=str2double(get(handles.checkbox_sSVD_th,'String'));
deconv_methods.cSVD.en=get(handles.checkbox_cSVD,'Value');
deconv_methods.cSVD.th=str2double(get(handles.checkbox_cSVD_th,'String'));
deconv_methods.oSVD.en=get(handles.checkbox_oSVD,'Value');
deconv_methods.oSVD.OI=str2double(get(handles.checkbox_oSVD_OI,'String'));
deconv_methods.tikhonov.en=get(handles.checkbox_tikhonov,'Value');
[CBF,CBV,MTT]=calc_maps(concentration4D,AIF,Mask,deconv_methods,handles);
maps_path=save_maps(CBF,CBV,MTT,deconv_methods);
maps_foler_message_str=['Maps saved to folder: ',maps_path];
set(handles.maps_folder_message,'Enable','on');
set(handles.maps_folder_message,'String',maps_foler_message_str);
set(handles.open_maps_folder_button,'Enable','On');
guidata(hObject,handles);

% --- Executes on button press in checkbox_cSVD.
function checkbox_cSVD_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_cSVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
cSVD_enable=get(hObject,'Value');
handles.deconv_methods.cSVD=cSVD_enable;
if cSVD_enable
   set(handles.checkbox_cSVD_th,'Enable','on');
   set(handles.text_cSVD_th,'Enable','on');
else
    set(handles.checkbox_cSVD_th,'Enable','off');
    set(handles.text_cSVD_th,'Enable','off');
end
guidata(hObject,handles);


% --- Executes on button press in checkbox_tikhonov.
function checkbox_tikhonov_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_tikhonov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.deconv_methods.tikhonov=get(hObject,'Value');
guidata(hObject,handles);

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function checkbox_sSVD_th_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_sSVD_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of checkbox_sSVD_th as text
%        str2double(get(hObject,'String')) returns contents of checkbox_sSVD_th as a double


% --- Executes during object creation, after setting all properties.
function checkbox_sSVD_th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_sSVD_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles=guidata(hObject);

guidata(hObject,handles);

% --- Executes on button press in checkbox_sSVD.
function checkbox_sSVD_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_sSVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
sSVD_enable=get(hObject,'Value');
handles.deconv_methods.sSVD=sSVD_enable;
if sSVD_enable
   set(handles.checkbox_sSVD_th,'Enable','on');
   set(handles.text_sSVD_th,'Enable','on');
else
    set(handles.checkbox_sSVD_th,'Enable','off');
    set(handles.text_sSVD_th,'Enable','off');
end
guidata(hObject,handles);


% --- Executes on button press in checkbox_oSVD.
function checkbox_oSVD_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_oSVD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
oSVD_enable=get(hObject,'Value');
handles.deconv_methods.oSVD=oSVD_enable;
if oSVD_enable
   set(handles.checkbox_oSVD_OI,'Enable','on');
   set(handles.text_oSVD_OI,'Enable','on');
else
    set(handles.checkbox_oSVD_OI,'Enable','off');
    set(handles.text_oSVD_OI,'Enable','off');
end
guidata(hObject,handles);



function checkbox_cSVD_th_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_cSVD_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of checkbox_cSVD_th as text
%        str2double(get(hObject,'String')) returns contents of checkbox_cSVD_th as a double


% --- Executes during object creation, after setting all properties.
function checkbox_cSVD_th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_cSVD_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function checkbox_oSVD_OI_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_oSVD_OI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of checkbox_oSVD_OI as text
%        str2double(get(hObject,'String')) returns contents of checkbox_oSVD_OI as a double


% --- Executes during object creation, after setting all properties.
function checkbox_oSVD_OI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_oSVD_OI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in deconvolution_status_textbox.
function deconvolution_status_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to deconvolution_status_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns deconvolution_status_textbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from deconvolution_status_textbox


% --- Executes during object creation, after setting all properties.
function deconvolution_status_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deconvolution_status_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_ZeroValues.
function checkbox_ZeroValues_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ZeroValues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ZeroValues


% --- Executes on button press in open_maps_folder_button.
function open_maps_folder_button_Callback(hObject, eventdata, handles)
% hObject    handle to open_maps_folder_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function maps_folder_message_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maps_folder_message (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
