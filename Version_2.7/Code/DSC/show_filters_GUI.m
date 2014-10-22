function varargout = show_filters_GUI(varargin)
% SHOW_FILTERS_GUI MATLAB code for show_filters_GUI.fig
%      SHOW_FILTERS_GUI, by itself, creates a new SHOW_FILTERS_GUI or raises the existing
%      singleton*.
%
%      H = SHOW_FILTERS_GUI returns the handle to a new SHOW_FILTERS_GUI or the handle to
%      the existing singleton*.
%
%      SHOW_FILTERS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOW_FILTERS_GUI.M with the given input arguments.
%
%      SHOW_FILTERS_GUI('Property','Value',...) creates a new SHOW_FILTERS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before show_filters_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to show_filters_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help show_filters_GUI

% Last Modified by GUIDE v2.5 27-Apr-2014 16:51:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @show_filters_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @show_filters_GUI_OutputFcn, ...
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


% --- Executes just before show_filters_GUI is made visible.
function show_filters_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to show_filters_GUI (see VARARGIN)

% Choose default command line output for show_filters_GUI
handles=guidata(hObject);
handles.output = hObject;

hMainGUI=getappdata(0,'hMainGUI');
data4D=getappdata(0,'data4D_init');
Mask=getappdata(0,'Mask');
BrainMask=getappdata(0,'BrainMask');
data_voxels_row_col_slice=getappdata(0,'data_voxels_row_col_slice');

num_brain_voxels=length(find(BrainMask>0));
num_data_voxels=size(data_voxels_row_col_slice,1);
set(handles.text1,'String',['View mask - ',num2str(round(100*num_data_voxels/num_brain_voxels*10)/10),'% of brain voxels passed.']);
%combine data4d with mask to one image where the data voxels are colored
%in red
show_data_time_sample=50;
red_vec=zeros(1,1,3);
red_vec(1)=1;
data_to_show=repmat(data4D(:,:,:,show_data_time_sample),[1,1,1,1,3]);
% data_to_show(Mask>0)=
handles.apply_filter=0;

axes(handles.show_filters_fig);
N_slices=size(data4D,3);
Nrows_subplot=ceil(sqrt(N_slices));
Ncols_subplot=ceil(N_slices/Nrows_subplot);
for slice=1:N_slices
    slice_data_to_show=squeeze(data_to_show(:,:,slice,:,:));
    max_val=max(max(slice_data_to_show(:,:,1)));
    slice_mask_inds=find(Mask(:,:,slice)>0);
    % marking the mask (on the slice) in red
    for color=1:3
        slice_color=slice_data_to_show(:,:,color);
        if color==1
            slice_color(slice_mask_inds)=max_val;
        else
            slice_color(slice_mask_inds)=0;
        end
        slice_data_to_show(:,:,color)=slice_color;
    end
    subplot(Nrows_subplot,Ncols_subplot,slice);
    imshow(mat2gray(slice_data_to_show));
    title(['Slice ',num2str(slice)]);
end


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes show_filters_GUI wait for user response (see UIRESUME)
% uiwait(handles.show_filters_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = show_filters_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in apply_filter_pushbutton.
function apply_filter_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to apply_filter_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.apply_filter=1;
guidata(hObject, handles);
setappdata(0,'apply_filter',1);
uiresume(gcbf);
% close(handles.show_filters_GUI)


% --- Executes on button press in back_pushbutton.
function back_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to back_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.apply_filter=0;
guidata(hObject, handles);
setappdata(0,'apply_filter',0);
uiresume(gcbf);
% close(handles.show_filters_GUI)
