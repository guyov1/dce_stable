function varargout = show_filters(varargin)
% SHOW_FILTERS MATLAB code for show_filters.fig
%      SHOW_FILTERS, by itself, creates a new SHOW_FILTERS or raises the existing
%      singleton*.
%
%      H = SHOW_FILTERS returns the handle to a new SHOW_FILTERS or the handle to
%      the existing singleton*.
%
%      SHOW_FILTERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOW_FILTERS.M with the given input arguments.
%
%      SHOW_FILTERS('Property','Value',...) creates a new SHOW_FILTERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before show_filters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to show_filters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help show_filters

% Last Modified by GUIDE v2.5 27-Apr-2014 16:14:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @show_filters_OpeningFcn, ...
                   'gui_OutputFcn',  @show_filters_OutputFcn, ...
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


% --- Executes just before show_filters is made visible.
function show_filters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to show_filters (see VARARGIN)

% Choose default command line output for show_filters
handles.output = hObject;

axes(handles.show_filters_fig);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes show_filters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = show_filters_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on button press in back_pushbutton.
function back_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to back_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
