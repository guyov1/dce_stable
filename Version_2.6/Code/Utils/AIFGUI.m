function varargout = AIFGUI(varargin)
% AIFGUI MATLAB code for AIFGUI.fig
%      AIFGUI, by itself, creates a new AIFGUI or raises the existing
%      singleton*.
%
%      H = AIFGUI returns the handle to a new AIFGUI or the handle to
%      the existing singleton*.
%
%      AIFGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AIFGUI.M with the given input arguments.
%
%      AIFGUI('Property','Value',...) creates a new AIFGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AIFGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AIFGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AIFGUI

% Last Modified by GUIDE v2.5 10-Nov-2013 10:57:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AIFGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @AIFGUI_OutputFcn, ...
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


% --- Executes just before AIFGUI is made visible.
function AIFGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AIFGUI (see VARARGIN)

% Choose default command line output for AIFGUI
handles.output = hObject;
% handles.BaseP='/home/john/Projects/code/john/Database/DCEOut/HaHa_20121212/';
% handles.BaseP='/media/OS/DATA/ArSi_20081130/';
% handles.BaseP='\\fmri-t9\users\Moran\lesionVasClassification\HealthyControls_New\GBPatients\DCEGUI_out\WhSa_20070813\';
handles.BaseP=varargin{1};

handles.Sliders=[handles.slider1 handles.slider2 handles.slider3 handles.slider4 handles.slider5 handles.slider6 handles.slider7 handles.slider8 handles.slider9];
handles.Titles=[handles.text1 handles.text2 handles.text3 handles.text4 handles.text5 handles.text6 handles.text7 handles.text8 handles.text9];
handles.Edits=[handles.edit1 handles.edit2 handles.edit3 handles.edit4 handles.edit5 handles.edit6 handles.edit7 handles.edit8 handles.edit9];
handles.Axes=[handles.axes2 handles.axes3 handles.axes4 handles.axes5 handles.axes6 handles.axes7 handles.axes8 handles.axes9 handles.axes10];
handles.nAxes=numel(handles.Axes);
Tmp=load([handles.BaseP 'PKM.mat']);
handles.DataToFit=Tmp.DataToFit2;
handles.GoodTs=Tmp.GoodTs;
handles.BaseAIFParams=Tmp.OutAIFParam;
handles.BaseTimeSamples = Tmp.SampleTs;
handles.InspectedAIFParamsFN=[handles.BaseP 'InspectedAIFParams.mat'];
handles.InspectedAIFParamsTimeFN=[handles.BaseP 'InspectedAIFParamsTime.mat'];
handles.CurAIFParams=Tmp.OutAIFParam;
if(exist(handles.InspectedAIFParamsFN,'file'))
    Tmpx=load(handles.InspectedAIFParamsFN);
    handles.CurAIFParams=Tmpx.InspectedParams;
end
handles.TimeBetweenDCEVolsFinal=Tmp.TimeBetweenDCEVolsFinal;
handles.nVols=Tmp.nSVols;
handles.MaxAmp=Tmp.MaxAmp;

InterpolationFactor=ceil(handles.TimeBetweenDCEVolsFinal);
DTimeVec=diff(handles.GoodTs);
F=find(DTimeVec>DTimeVec(1)*2,1);
handles.nNormalTs=numel(handles.GoodTs);
if(~isempty(F))
    handles.nNormalTs=F(1);
end

HInterpolationFactor=ceil(InterpolationFactor*Tmp.Options.SubSecondResolution);
Hdt=handles.TimeBetweenDCEVolsFinal/(60*HInterpolationFactor);
handles.HSampleTs=[0:Hdt:handles.GoodTs(handles.nNormalTs) handles.GoodTs((handles.nNormalTs+1):end)];
handles.nHVols=numel(handles.HSampleTs);
handles.Hdt=Hdt;

ThreeSec=ceil(Tmp.Options.MaxTDif_ForAIFSearch/(Hdt*60));
handles.TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;
handles.nTDif=numel(handles.TDif);

% handles.MaxAmp=max(handles.DataToFit(:));
handles.AIF_Parker9t=@(x,t) AIF_Parkerg3( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8),x(9),max(handles.GoodTs))*x(2);
% handles.AIFFunc=@(x,t) handles.AIF_Parker9t([x(1:8) x(9)/handles.MaxAmp],t).*handles.MaxAmp;
handles.AIFFunc=@(x,t) handles.AIF_Parker9t(x,t);
handles.TitleStr={'T1','Amp1','Std1','T2','Amp2','Std2','SlopeHeight','SlopeExp','Base level'};
handles.CurIdxs=randsample(size(handles.DataToFit,1),handles.nAxes,true);
for i=1:numel(handles.Titles)
    set(handles.Titles(i),'String',handles.TitleStr{i});
end
MinFirstBolusSig=1;
LB=[0.1 0     MinFirstBolusSig/60 0.1     0   0.01  0.0001   0  0]';
UB=[10  1.5   0.25                10      1   1     2        2  1]';
handles.LB=LB;
handles.UB=UB;
for i=1:numel(handles.Sliders)
    set(handles.Sliders(i),'Min',handles.LB(i));
    set(handles.Sliders(i),'Max',handles.UB(i));
%     set(handles.Sliders(i),'SliderStep',[(handles.UB(i)-handles.LB(i))/10 (handles.UB(i)-handles.LB(i))/100]);
    set(handles.Sliders(i),'SliderStep',[1/100 1/10]);
    set(handles.Sliders(i),'Value',handles.CurAIFParams(i));
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AIFGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AIFGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
HandleSliders(hObject, eventdata, handles);

function ParamChange(hObject, eventdata, handles)
handles.HAIF=handles.AIFFunc(handles.CurAIFParams,handles.HSampleTs);
CHAIF=cumtrapz(handles.HSampleTs,handles.HAIF);
CHAIF=[0 cumsum(handles.HAIF(1:end-1))]*handles.Hdt;

SAIF=zeros([handles.nTDif numel(handles.GoodTs)]);
CSAIF=zeros([handles.nTDif numel(handles.GoodTs)]);
for i=1:handles.nTDif
    SAIF(i,:)=interp1(handles.HSampleTs,handles.HAIF,handles.GoodTs+handles.TDif(i),[],'extrap');
    CSAIF(i,:)=interp1(handles.HSampleTs,CHAIF,handles.GoodTs+handles.TDif(i),[],'extrap');
end
[handles.PKs handles.Sims] = FindPKBATgAIFMuraseF(handles.DataToFit,SAIF,handles.GoodTs,CSAIF);

guidata(hObject, handles);
ShowIm(hObject, eventdata, handles);

function ShowIm(hObject, eventdata, handles)
for i=1:numel(handles.Sliders)
    set(handles.Sliders(i),'Value',handles.CurAIFParams(i));
end
axes(handles.axes1);cla;
plot(handles.HSampleTs,handles.HAIF,'r-');
hold on;
plot(handles.GoodTs,handles.AIFFunc(handles.CurAIFParams,handles.GoodTs),'b*');

for i=1:numel(handles.Axes)
    axes(handles.Axes(i));cla;
    plot(handles.GoodTs,handles.DataToFit(handles.CurIdxs(i),:),'k*');
    hold on;
    plot(handles.GoodTs,handles.Sims(handles.CurIdxs(i),:),'b-');
end

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
HandleSliders(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
HandleSliders(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
HandleSliders(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
HandleSliders(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
HandleSliders(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
HandleSliders(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
HandleSliders(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider8_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
HandleSliders(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
HandleEdits(hObject, eventdata, handles);

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



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
HandleEdits(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
HandleEdits(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
HandleEdits(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
HandleEdits(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
HandleEdits(hObject, eventdata, handles);

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



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
HandleEdits(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
HandleEdits(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function HandleEdits(hObject, eventdata, handles)
for i=1:numel(handles.Sliders)
    CurV=str2num(get(handles.Edits(i),'String'));
    set(handles.Sliders(i),'Value',CurV);
    handles.CurAIFParams(i)=CurV;
end
guidata(hObject, handles);
ParamChange(hObject, eventdata, handles);

function HandleSliders(hObject, eventdata, handles)
for i=1:numel(handles.Sliders)
    CurV=get(handles.Sliders(i),'Value');
    set(handles.Edits(i),'String',num2str(CurV));
    handles.CurAIFParams(i)=CurV;
end
guidata(hObject, handles);
ParamChange(hObject, eventdata, handles);



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
HandleEdits(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider9_Callback(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
HandleSliders(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function slider9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushRevert.
function pushRevert_Callback(hObject, eventdata, handles)
% hObject    handle to pushRevert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.CurAIFParams   =handles.BaseAIFParams;
handles.CurTimeSamples =handles.BaseTimeSamples;
for i=1:numel(handles.Sliders)
    set(handles.Edits(i),'String',num2str(handles.BaseAIFParams(i)));
end
guidata(hObject, handles);
HandleEdits(hObject, eventdata, handles);

% --- Executes on button press in pushSave.
function pushSave_Callback(hObject, eventdata, handles)
% hObject    handle to pushSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
InspectedParams = handles.CurAIFParams;
InspectedParamsTimeSamples = handles.BaseTimeSamples;

save(handles.InspectedAIFParamsFN,'InspectedParams');
save(handles.InspectedAIFParamsTimeFN,'InspectedParamsTimeSamples');



disp('Inspected AIF parameters saved.');