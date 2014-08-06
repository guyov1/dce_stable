function varargout = DCEResultsGUI(varargin)
% DCERESULTSGUI M-file for DCEResultsGUI.fig
%      DCERESULTSGUI, by itself, creates a new DCERESULTSGUI or raises the existing
%      singleton*.
%
%      H = DCERESULTSGUI returns the handle to a new DCERESULTSGUI or the handle to
%      the existing singleton*.
%
%      DCERESULTSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DCERESULTSGUI.M with the given input arguments.
%
%      DCERESULTSGUI('Property','Value',...) creates a new DCERESULTSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DCEResultsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DCEResultsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DCEResultsGUI

% Last Modified by GUIDE v2.5 25-Jul-2012 12:22:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DCEResultsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DCEResultsGUI_OutputFcn, ...
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


% --- Executes just before DCEResultsGUI is made visible.
function DCEResultsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DCEResultsGUI (see VARARGIN)

% Choose default command line output for DCEResultsGUI
handles.output = hObject;

handles.Vols=varargin{1};
handles.Titles=varargin{2};
handles.Msk=varargin{3};
handles.CTC2D=varargin{4};
handles.Idx3D=varargin{5};
% handles.MIdxs=varargin{6};
handles.PKs=varargin{7};
handles.HSAIF=varargin{8};
handles.HSHConvd=varargin{9};
handles.Keps1I=varargin{10};
handles.SampleTs=varargin{11};
handles.HSampleTs=varargin{12};
handles.HAIF=varargin{13};

handles.nSlices=size(handles.Vols{1},3);
set(handles.slider1,'Min',1);
set(handles.slider1,'Max',handles.nSlices);
set(handles.slider1,'Value',floor(handles.nSlices/2));
set(handles.slider1,'SliderStep',ones(1,2)*(1/(handles.nSlices-1)));

set(handles.listbox1,'String',handles.Titles);
set(handles.listbox1,'Value',1);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DCEResultsGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DCEResultsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
ShowIm(hObject, eventdata, handles);

axes(handles.axes2);cla;
plot(handles.HSampleTs,handles.HAIF,'k');

function ShowIm(hObject, eventdata, handles)
CurV=get(handles.listbox1,'Value');
CurS=round(get(handles.slider1,'Value'));

axes(handles.axes1);

Data=squeeze(handles.Vols{CurV}(:,:,CurS));
% [out minval maxval]=normalP(Data,0, 0.95);
if(numel(unique(Data))>2)
    [out minval maxval]=normalP(handles.Vols{CurV},0.03, 0.95);
else
    minval=0;
    maxval=1;
end
imagesc(mritransform((Data)),[minval maxval]);
colorbar;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
% title([num2str(

% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
ShowIm(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
ShowIm(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ChosenVoxels=ginput;
SDCE=size(handles.Vols{1});
CurS=get(handles.slider1,'Value');
ChosenVoxels=ceil([ChosenVoxels(:,1) SDCE(2)+1-ChosenVoxels(:,2) ChosenVoxels(:,1)*0+CurS]);
Idxs=handles.Idx3D(sub2ind(SDCE,ChosenVoxels(:,1),ChosenVoxels(:,2),ChosenVoxels(:,3)));
Good=isfinite(Idxs);
Idxs=Idxs(Good);
ChosenVoxels=ChosenVoxels(Good,:);

axes(handles.axes2);cla;

if(isempty(Idxs))
    plot(handles.HSampleTs,handles.HAIF,'k');
    set(handles.listbox1,'String',handles.Titles);
    return;
end

plot(handles.SampleTs,handles.CTC2D(Idxs,:)','*');
Strs=handles.Titles;

for i=1:numel(Idxs)
    Regressors=[handles.HSAIF(handles.MIdxs(Idxs(i),1),:); squeeze(handles.HSHConvd(handles.MIdxs(Idxs(i),1),handles.Keps1I(handles.MIdxs(Idxs(i),2)),:))'];
    Sims(i,:)=((Regressors')*handles.CXs(:,Idxs(i)));
    for j=1:numel(handles.Titles)
        Strs{j}=[Strs{j} ' ' num2str(handles.Vols{j}(ChosenVoxels(i,1),ChosenVoxels(i,2),ChosenVoxels(i,3)))];
    end
end
set(handles.listbox1,'String',Strs);
hold on;
plot(handles.HSampleTs,Sims','--','LineWidth',0.5);
title([num2str(ChosenVoxels(1,1)) ' ' num2str(ChosenVoxels(1,2)) ' ' num2str(CurS)]);