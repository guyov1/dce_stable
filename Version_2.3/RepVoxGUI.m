function varargout = RepVoxGUI(varargin)
% REPVOXGUI MATLAB code for RepVoxGUI.fig
%      REPVOXGUI, by itself, creates a new REPVOXGUI or raises the existing
%      singleton*.
%
%      H = REPVOXGUI returns the handle to a new REPVOXGUI or the handle to
%      the existing singleton*.
%
%      REPVOXGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REPVOXGUI.M with the given input arguments.
%
%      REPVOXGUI('Property','Value',...) creates a new REPVOXGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RepVoxGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RepVoxGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RepVoxGUI

% Last Modified by GUIDE v2.5 22-Dec-2013 15:16:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RepVoxGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @RepVoxGUI_OutputFcn, ...
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


% --- Executes just before RepVoxGUI is made visible.
function RepVoxGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RepVoxGUI (see VARARGIN)

% Choose default command line output for RepVoxGUI
handles.output = hObject;

% handles.BaseP='/media/OS/STRDCE/John/Database/DCEOut/HaMo_20130714/';
handles.BaseP=varargin{1};
% handles.Vol=loadniidata([handles.BaseP 'Baseline.nii']);
D=dir([handles.BaseP 'DCEMainCoreged' filesep '*.nii']);
handles.Vol=loadniidata([handles.BaseP 'DCEMainCoreged' filesep D(end).name]);
InspectedRepVoxNII=[handles.BaseP 'InspectedRepVox.nii'];
if(exist(InspectedRepVoxNII,'file'))
    Tmp=loadniidata([handles.BaseP 'InspectedRepVox.nii']);
else
    Tmp=loadniidata([handles.BaseP 'ChosenVoxelsForAIFFinding.nii']);
end
handles.Sz=size(handles.Vol);
F=find(Tmp);
[I{1} I{2} I{3}]=ind2sub(handles.Sz,F);
handles.I=[I{1} I{2} I{3}];
Tmp=num2str(handles.I);
Tmp=mat2cell(Tmp,ones(size(Tmp,1),1),size(Tmp,2));
set(handles.listbox1,'String',Tmp);
set(handles.listbox1,'Value',[]);

handles.V4D=loadniidata([handles.BaseP 'CTC4D.nii']);

[Tmp handles.DR]=FindDR(handles.Vol);

handles.nSlices=size(handles.Vol,3);
set(handles.slider1,'Min',1);
set(handles.slider1,'Max',handles.nSlices);
set(handles.slider1,'Value',floor(handles.nSlices/2));
set(handles.slider1,'SliderStep',ones(1,2)*(1/(handles.nSlices-1)));

WorkingP=handles.BaseP;
PrepareFN=[WorkingP 'AfterPrepare4D.mat'];
load(PrepareFN,'Options');
UnderSampling=Options.SubSampling;
if(UnderSampling==1)
    USStr='';
else
    USStr=['_' num2str(UnderSampling)];
end
CTCFN=[WorkingP 'AfterCTC' USStr '.mat'];
load(CTCFN,'NumVols','TimeBetweenDCEVolsFinal','Options','BolusStart')
handles.Options=Options;
handles.NumVols=NumVols;
handles.TimeBetweenDCEVolsFinal=TimeBetweenDCEVolsFinal;
handles.BolusStart=BolusStart;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RepVoxGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RepVoxGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
ShowIm(hObject, eventdata, handles);

function ShowIm(hObject, eventdata, handles)
CurV=get(handles.listbox1,'Value');
CurS=round(get(handles.slider1,'Value'));

axes(handles.axes1);

RGB=repmat(handles.Vol(:,:,CurS),[1 1 3]);
RGB=RGB/handles.DR;
RGB=min(RGB,1);
CurVox=handles.I(handles.I(:,3)==CurS,[1 2]);
Msk=squeeze(handles.Vol(:,:,CurS))*0>1;
Msk2=Msk;
for i=1:size(CurVox,1)
    Msk(CurVox(i,1),CurVox(i,2))=true;
end
for i=1:numel(CurV)
    if(handles.I(CurV(i),3)==CurS)
        Msk2(handles.I(CurV(i),1),handles.I(CurV(i),2))=true;
    end
end
Msk3=false(size(handles.Vol));
if(isfield(handles,'ChosenVoxel'))
    Msk3(handles.ChosenVoxel(1),handles.ChosenVoxel(2),handles.ChosenVoxel(3))=true;
end
Msk3=squeeze(Msk3(:,:,CurS));
TmpD=imdilate(Msk,strel('disk',2));
TmpD2=imdilate(Msk2,strel('disk',2));
TmpD3=imdilate(Msk3,strel('disk',2));
TmpD4=false(size(handles.Vol));
if(isfield(handles,'ROI'))
    TmpD4=handles.ROI;
end
TmpD4=squeeze(TmpD4(:,:,CurS));
for i=1:3
    Tmp=squeeze(RGB(:,:,i));
    if((i==1) || (i==3))
        Alpha=0.5;
        Tmp(TmpD4)=Alpha+(1-Alpha)*Tmp(TmpD4);
    end
    Tmp(TmpD)=(i==3);
    Tmp(TmpD3)=(i==1) || (i==2);
    Tmp(Msk)=(i==1);
    Tmp(TmpD2)=(i==2);
    Tmp(Msk2)=(i==1);
    Tmp(Msk3)=(i==1);
    RGB(:,:,i)=Tmp;
end

imagesc(mritransform(RGB));
title(CurS);

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

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
CurV=get(handles.listbox1,'Value');
ChosenVoxel=handles.I(CurV,:);
CurS=handles.I(CurV(end),3);
set(handles.slider1,'Value',CurS);
axes(handles.axes2);cla;
Tmp=[];
for i=1:size(ChosenVoxel,1)
    Tmp(:,i)=squeeze(handles.V4D(ChosenVoxel(i,1),ChosenVoxel(i,2),ChosenVoxel(i,3),:));
end
plot(Tmp);
guidata(hObject, handles);
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

% --- Executes on button press in pushChoose.
function pushChoose_Callback(hObject, eventdata, handles)
% hObject    handle to pushChoose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ChosenVoxel=ginput(1);
CurS=get(handles.slider1,'Value');
ChosenVoxel=ceil([ChosenVoxel(1) handles.Sz(2)+1-ChosenVoxel(2)]);
handles.ChosenVoxel=[ChosenVoxel CurS];
CurVoxI=find(handles.I(:,3)==CurS);
CurVox=handles.I(CurVoxI,[1 2]);
[a b]=min(sum((CurVox-repmat(ChosenVoxel(:,1:2),[size(CurVox,1) 1])).^2,2));
set(handles.listbox1,'Value',CurVoxI(b));
axes(handles.axes2);cla;
plot(squeeze(handles.V4D(ChosenVoxel(1),ChosenVoxel(2),CurS,:)));
guidata(hObject, handles);
ShowIm(hObject, eventdata, handles);

% --- Executes on button press in pushAdd.
function pushAdd_Callback(hObject, eventdata, handles)
% hObject    handle to pushAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~isfield(handles,'ChosenVoxel'))
    return;
end
ChosenVoxel=handles.ChosenVoxel;
% CurS=get(handles.slider1,'Value');
CurS=ChosenVoxel(3);
IA=handles.I(:,3)<=CurS;
IB=handles.I(:,3)>CurS;
% handles.I=[handles.I(IA,:); [ChosenVoxel CurS]; handles.I(IB,:)];
handles.I=[handles.I(IA,:); ChosenVoxel; handles.I(IB,:)];
Tmp=num2str(handles.I);
Tmp=mat2cell(Tmp,ones(size(Tmp,1),1),size(Tmp,2));
set(handles.listbox1,'String',Tmp);
set(handles.listbox1,'Value',sum(IA)+1);
guidata(hObject, handles);

ShowIm(hObject, eventdata, handles);

% --- Executes on button press in pushRemove.
function pushRemove_Callback(hObject, eventdata, handles)
% hObject    handle to pushRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
V=get(handles.listbox1,'Value');
handles.I(V,:)=[];
Tmp=num2str(handles.I);
Tmp=mat2cell(Tmp,ones(size(Tmp,1),1),size(Tmp,2));
set(handles.listbox1,'String',Tmp);
set(handles.listbox1,'Value',[]);

if(isfield(handles,'ChosenVoxel'))
    ChosenVoxel=handles.ChosenVoxel;
%     CurS=get(handles.slider1,'Value');
    CurS=ChosenVoxel(3);
    CurVoxI=find(handles.I(:,3)==CurS);
    CurVox=handles.I(CurVoxI,[1 2]);
    [a b]=min(sum((CurVox-repmat(ChosenVoxel(1:2),[size(CurVox,1) 1])).^2,2));
    set(handles.listbox1,'Value',CurVoxI(b));
end
guidata(hObject, handles);

ShowIm(hObject, eventdata, handles);


% --- Executes on button press in pushSave.
function pushSave_Callback(hObject, eventdata, handles)
% hObject    handle to pushSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Tmp=handles.Vol*0>1;
for i=1:size(handles.I)
    Tmp(handles.I(i,1),handles.I(i,2),handles.I(i,3))=true;
end
Raw2Nii(Tmp*1,[handles.BaseP 'InspectedRepVox.nii'],'float32', [handles.BaseP 'DCEMean.nii']);
disp(['Inspected representative voxels saved to ' handles.BaseP 'InspectedRepVox.nii']);


% --- Executes on button press in pushROI.
function pushROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurS=get(handles.slider1,'Value');
Ps=gginput;
ChosenVoxels=ceil([Ps(:,1) handles.Sz(2)+1-Ps(:,2)]);
I=squeeze(handles.Vol(:,:,CurS));
BW = roipoly(I,ChosenVoxels(:,2),ChosenVoxels(:,1));
% figure;imagesc(mritransform(I+BW*100));
Msk=handles.Vol*0>1;
Msk(:,:,CurS)=BW>0;
handles.ROI=Msk;
CurData=Reshape4d22d(handles.V4D,Msk);
MskF=loadniidata([handles.BaseP 'ArtMsk.nii']);
MskOut=MskF & Msk;
% [DataOut, Idxs, MskOut]=FilterCTCToFindArt(CurData,Msk,handles.NumVols,handles.TimeBetweenDCEVolsFinal,handles.BolusStart,handles.Options);
[I J K]=ind2sub(size(MskOut),find(MskOut));
ChosenVoxels=[I J];
handles.I=[I J K];
Tmp=num2str(handles.I);
Tmp=mat2cell(Tmp,ones(size(Tmp,1),1),size(Tmp,2));
set(handles.listbox1,'String',Tmp);
set(handles.listbox1,'Value',[]);
set(handles.listbox1,'ListboxTop',1);
guidata(hObject, handles);

ShowIm(hObject, eventdata, handles);
