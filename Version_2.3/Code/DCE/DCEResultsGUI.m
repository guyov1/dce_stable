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

% Last Modified by GUIDE v2.5 10-Nov-2013 11:30:54

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
handles.TTL=varargin{14};
handles.TDif=varargin{15};
handles.nTDif=numel(varargin{15});

CMaps.MeanVol='cool';
CMaps.BAT='gray';
try
    Tmp=getComputerParams('CMaps');
    CMaps=Tmp;
catch
    setComputerParamM('CMaps',CMaps)
end

handles.CMaps=CMaps;

for CurV=1:numel(handles.Vols)
    % [out minval maxval]=normalP(Data,0, 0.95);
%     if(numel(unique(handles.Vols{CurV}(isfinite(handles.Vols{CurV}))))>2)
%         %     [out minval maxval]=normalP(handles.Vols{CurV},0.03, 0.95);
%         [minval maxval]=FindDR(handles.Vols{CurV});
%     else
%         minval=0;
%         maxval=1;
%     end
    handles.DR(CurV,:)=[NaN NaN];
end

handles.nSlices=size(handles.Vols{1},3);
set(handles.slider1,'Min',1);
set(handles.slider1,'Max',handles.nSlices);
set(handles.slider1,'Value',floor(handles.nSlices/2));
set(handles.slider1,'SliderStep',ones(1,2)*(1/(handles.nSlices-1)));

set(handles.listbox1,'String',handles.Titles);
set(handles.listbox1,'Value',1);

set(handles.edit5,'String','0');
set(handles.edit6,'String',num2str(handles.SampleTs(end)));

LUTs=load('CronMaps.mat');
LUTs=LUTs.CronMaps;
handles.LUTs=LUTs;
set(handles.popupcmap,'String',cellfun(@(x) x(2:end),fieldnames(LUTs),'UniformOutput',false));

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
if(isnan(handles.DR(CurV,1)))
    if(numel(unique(handles.Vols{CurV}(isfinite(handles.Vols{CurV}))))>20)
        %     [out minval maxval]=normalP(handles.Vols{CurV},0.03, 0.95);
        [minval maxval]=FindDR(handles.Vols{CurV});
    else
        minval=min(handles.Vols{CurV}(isfinite(handles.Vols{CurV})));
        maxval=max(handles.Vols{CurV}(isfinite(handles.Vols{CurV})));
    end
    if(isempty(minval))
        handles.DR(CurV,:)=[0 1];
    else
        handles.DR(CurV,:)=[minval maxval];
    end
    guidata(hObject, handles);
end
set(handles.edit1,'String',num2str(handles.DR(CurV,1)));
set(handles.edit2,'String',num2str(handles.DR(CurV,2)));
%Data(isnan(Data))=NaN;
imagesc(mritransform(Data),handles.DR(CurV,:));
set(gca,'XTick',[]);
set(gca,'YTick',[]);
title(strrep(handles.TTL,'\','\\'));

UseRGB=false;
if(isfield(handles,'ChosenVoxels'))
    for i=1:size(handles.ChosenVoxels,1)
        if(handles.ChosenVoxels(i,3)==CurS)
            UseRGB=true;
        end
    end
end
% if(UseRGB)
%     colormap(gray);
%
ImagesS=get(handles.listbox1,'String');
% CurIName=ImagesS{CurV};
CurIName=handles.Titles{CurV}; %Tmp{CurV};
CMapsS=get(handles.popupcmap,'String');
CMaps=handles.CMaps;
if(isfield(CMaps,CurIName))
    colormap(handles.LUTs.(['A' CMaps.(CurIName)]));
    set(handles.popupcmap,'Value',find(strcmp(CMapsS,CMaps.(CurIName))));
else
    colormap('default');
    set(handles.popupcmap,'Value',1);
end
% end
CurCMAP=jet(256);
if(isfield(CMaps,CurIName))
    CurCMAP=handles.LUTs.(['A' CMaps.(CurIName)]);
end
SliceData=min(1,max(0,(mritransform(Data)-handles.DR(CurV,1))./diff(handles.DR(CurV,:))));
RGB=ind2rgb(floor(SliceData*256),CurCMAP);
Datab=squeeze(handles.Vols{1}(:,:,CurS));
SliceDatab=min(1,max(0,(mritransform(Datab)-handles.DR(1,1))./diff(handles.DR(1,:))));
RGBb=ind2rgb(floor(SliceDatab*256),gray(256));
TmpMsk=SliceData==0;
for i=1:3
    Tmp=squeeze(RGB(:,:,i));
    Tmpb=squeeze(RGBb(:,:,i));
    Tmp(TmpMsk)=Tmpb(TmpMsk);
    RGB(:,:,i)=Tmp;
end
imshow(RGB);
colorbar;

CLRs='kgrcmywwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww';
SDCE=size(handles.Vols{1});
if(isfield(handles,'ChosenVoxels'))
    for i=1:size(handles.ChosenVoxels,1)
        if(handles.ChosenVoxels(i,3)==CurS)
            rectangle('Position',[handles.ChosenVoxels(i,1)-5, SDCE(2)+1-handles.ChosenVoxels(i,2)-5,10,10],'LineWidth',2,'EdgeColor',CLRs(i));
        end
    end
end

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
ChosenVoxels=ChosenVoxels(all(ChosenVoxels>0,2),:);
SDCE=size(handles.Vols{1});
CurS=get(handles.slider1,'Value');
ChosenVoxels=ceil([ChosenVoxels(:,1) SDCE(2)+1-ChosenVoxels(:,2) ChosenVoxels(:,1)*0+CurS]);
Idxs=handles.Idx3D(sub2ind(SDCE,ChosenVoxels(:,1),ChosenVoxels(:,2),ChosenVoxels(:,3)));
Good=isfinite(Idxs);
Idxs=Idxs(Good);
ChosenVoxels=ChosenVoxels(Good,:);
handles.ChosenVoxels=ChosenVoxels;
guidata(hObject, handles);
PlotStuff(hObject, eventdata, handles);

function PlotStuff(hObject, eventdata, handles)
SDCE=size(handles.Vols{1});
ChosenVoxels=handles.ChosenVoxels;
Idxs=handles.Idx3D(sub2ind(SDCE,ChosenVoxels(:,1),ChosenVoxels(:,2),ChosenVoxels(:,3)));
Good=isfinite(Idxs);
Idxs=Idxs(Good);
CurS=get(handles.slider1,'Value');

axes(handles.axes2);cla;

if(isempty(Idxs))
    plot(handles.HSampleTs,handles.HAIF,'k');
    set(handles.listbox1,'String',handles.Titles);
    ax=axis();
    axis([ax(1:2) min(handles.HAIF) max(handles.HAIF)]);
    handles.ChosenVoxels=zeros(0,3);
    guidata(hObject, handles);
    ShowIm(hObject, eventdata, handles);
    return;
end

NormalizeCTCs=get(handles.checkNormalizeRows(1),'Value');

% NormalizeCTCs=1;

Strs=handles.Titles;

HConvIdxM=CreateConvIdxMFromSampleTs(numel(handles.HSampleTs));
HTriB=HConvIdxM>0;
HConvIdxMTriB=HConvIdxM(HTriB);

% BATfinal VpFinal KtransFinal Kepfinal VeFinal
BATIdx=1;
VpIdx=2;
KtransIdx=3;
KepIdx=4;
CurKeps=handles.PKs(Idxs,KepIdx);
CurKeps(isnan(CurKeps))=0;

Hdt=diff(handles.HSampleTs);
Hdt=Hdt(1);

HHSampleTs=0:Hdt:handles.HSampleTs(end);
HHAIF=interp1(handles.HSampleTs,handles.HAIF',HHSampleTs);
HHConvIdxM=CreateConvIdxMFromSampleTs(numel(HHSampleTs));
HHTriB=HHConvIdxM>0;
HHConvIdxMTriB=HHConvIdxM(HHTriB);

HHConvd2=DCECostFuncgrT1ForConv(HHAIF',CurKeps,HHSampleTs,HHConvIdxMTriB,HHTriB);
for i=1:size(HHConvd2,1)
    HConvd2(i,:)=interp1(HHSampleTs,HHConvd2(i,:),handles.HSampleTs,[],'extrap');
end
% HConvd2=DCECostFuncgrT1ForConv(handles.HAIF',CurKeps,handles.HSampleTs,HConvIdxMTriB,HTriB);

Hdt=handles.HSampleTs(2)-handles.HSampleTs(1);
dt=diff(handles.SampleTs(1:2));
% ThreeSec=ceil(3/(Hdt*60));
% TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;

CurBATs=handles.PKs(Idxs,BATIdx)/-60;
CurBATs(isnan(CurBATs))=1;
for i=1:numel(Idxs)
%     SHConvd2(i,:)=interp1(handles.HSampleTs,HConvd2(i,:)',handles.HSampleTs+handles.TDif(CurBATs(i)),[],'extrap')';
%     SHSAIF(i,:)=interp1(handles.HSampleTs,handles.HAIF,handles.HSampleTs+handles.TDif(CurBATs(i)),[],'extrap');
    SHConvd2(i,:)=interp1(handles.HSampleTs,HConvd2(i,:)',handles.HSampleTs+CurBATs(i),[],'extrap')';
    SHSAIF(i,:)=interp1(handles.HSampleTs,handles.HAIF,handles.HSampleTs+CurBATs(i),[],'extrap');
end
%
CurKtranses=handles.PKs(Idxs,KtransIdx);
CurKtranses(isnan(CurKtranses))=0;
for i=1:numel(Idxs)
%     Regressors=[SAIF(PKs(c,1),:); squeeze(SHConvd2(PKs(c,1),c,:))'];
    Regressors=[SHSAIF(i,:); squeeze(SHConvd2(i,:))];
        % BATfinal VpFinal KtransFinal Kepfinal VeFinal
    Sims(i,:)=((Regressors')*([handles.PKs(Idxs(i),[VpIdx]) CurKtranses(i)]'));        
%     figure;plot(handles.SampleTs,handles.CTC2D(Idxs,:),'k*',handles.SampleTs,Sims,'r',handles.SampleTs,SHSAIF*handles.PKs(Idxs(i),3),'g',handles.SampleTs,SHConvd2*handles.PKs(Idxs(i),4)/dt,'m')

%     Regressors=[handles.HSAIF(handles.MIdxs(Idxs(i),1),:); squeeze(handles.HSHConvd(handles.MIdxs(Idxs(i),1),handles.Keps1I(handles.MIdxs(Idxs(i),2)),:))'];
%     Sims(i,:)=((Regressors')*handles.CXs(:,Idxs(i)));
    for j=1:numel(handles.Titles)
        Strs{j}=[Strs{j} ' ' num2str(handles.Vols{j}(ChosenVoxels(i,1),ChosenVoxels(i,2),ChosenVoxels(i,3)))];
    end
end
set(handles.listbox1,'String',Strs);
if(NormalizeCTCs)
    plot(handles.HSampleTs,handles.HAIF./max(handles.HAIF),'k','LineWidth',3);hold on;
    Coeffs=max(Sims,[],2);
    Tmp=repMulti(handles.CTC2D(Idxs,:),1./Coeffs);
    plot(handles.SampleTs,Tmp','x');
else
    plot(handles.SampleTs,handles.CTC2D(Idxs,:)','x');
end
hold on;
if(NormalizeCTCs)
    plot(handles.HSampleTs,NormalizeByRows(Sims)','-','LineWidth',1);
else
    plot(handles.HSampleTs,Sims','-','LineWidth',1);
end
title([num2str(ChosenVoxels(1,1)) ' ' num2str(ChosenVoxels(1,2)) ' ' num2str(CurS)]);
ax=axis();
if(NormalizeCTCs)
    M=1;
else
    M=max(max(handles.CTC2D(Idxs,:)));
end
StartT=str2num(get(handles.edit5,'String'));
EndT=str2num(get(handles.edit6,'String'));
axis([StartT EndT -M/2 M*1.5]); %[-1 2]/1000


ShowIm(hObject, eventdata, handles);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
V=str2double(get(hObject,'String'));
if(~isfinite(V))
    return;
end
CurV=get(handles.listbox1,'Value');
if(V>=handles.DR(CurV,2))
    return;
end
handles.DR(CurV,1)=V;
guidata(hObject, handles);
ShowIm(hObject, eventdata, handles);


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
V=str2double(get(hObject,'String'));
if(~isfinite(V))
    return;
end
CurV=get(handles.listbox1,'Value');
if(V<=handles.DR(CurV,1))
    return;
end
handles.DR(CurV,2)=V;
guidata(hObject, handles);
ShowIm(hObject, eventdata, handles);


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


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurV=get(handles.listbox1,'Value');

% [out minval maxval]=normalP(Data,0, 0.95);
if(numel(unique(handles.Vols{CurV}(isfinite(handles.Vols{CurV}))))>2)
    %     [out minval maxval]=normalP(handles.Vols{CurV},0.03, 0.95);
    [minval maxval]=FindDR(handles.Vols{CurV});
else
    minval=0;
    maxval=1;
end
handles.DR(CurV,:)=[minval maxval];
guidata(hObject, handles);
ShowIm(hObject, eventdata, handles);


% --- Executes on button press in checkNormalizeRows.
function checkNormalizeRows_Callback(hObject, eventdata, handles)
% hObject    handle to checkNormalizeRows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkNormalizeRows
PlotStuff(hObject, eventdata, handles);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ChosenVoxels=ginput;
ChosenVoxels=ChosenVoxels(all(ChosenVoxels>0,2),:);
SDCE=size(handles.Vols{1});
CurS=get(handles.slider1,'Value');
ChosenVoxels=ceil([ChosenVoxels(:,1) SDCE(2)+1-ChosenVoxels(:,2) ChosenVoxels(:,1)*0+CurS]);
Idxs=handles.Idx3D(sub2ind(SDCE,ChosenVoxels(:,1),ChosenVoxels(:,2),ChosenVoxels(:,3)));
Good=isfinite(Idxs);
Idxs=Idxs(Good);
ChosenVoxels=ChosenVoxels(Good,:);
if(isfield(handles,'ChosenVoxels'))
    handles.ChosenVoxels=[handles.ChosenVoxels; ChosenVoxels];
else
    handles.ChosenVoxels=ChosenVoxels;
end
guidata(hObject, handles);
PlotStuff(hObject, eventdata, handles);


% --- Executes on selection change in popupcmap.
function popupcmap_Callback(hObject, eventdata, handles)
% hObject    handle to popupcmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CurM=get(handles.popupcmap,'Value');
CMapsS=get(handles.popupcmap,'String');
CurCMName=CMapsS{CurM};

CurV=get(handles.listbox1,'Value');
Tmp=get(handles.listbox1,'String');
CurIName=handles.Titles{CurV}; %Tmp{CurV};
CMaps=handles.CMaps;
CMaps.(CurIName)=CurCMName;
setComputerParamM('CMaps',CMaps);
handles.CMaps=CMaps;
guidata(hObject, handles);
ShowIm(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function popupcmap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupcmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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
PlotStuff(hObject, eventdata, handles);

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
PlotStuff(hObject, eventdata, handles);

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
