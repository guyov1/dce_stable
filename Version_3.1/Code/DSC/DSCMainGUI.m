function varargout = DSCMainGUI(varargin)
% DSCMAINGUI MATLAB code for DSCMainGUI.fig
%      DSCMAINGUI, by itself, creates a new DSCMAINGUI or raises the existing
%      singleton*.
%
%      H = DSCMAINGUI returns the handle to a new DSCMAINGUI or the handle to
%      the existing singleton*.
%
%      DSCMAINGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSCMAINGUI.M with the given input arguments.
%
%      DSCMAINGUI('Property','Value',...) creates a new DSCMAINGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DSCMainGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DSCMainGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DSCMainGUI

% Last Modified by GUIDE v2.5 09-Sep-2014 14:20:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DSCMainGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DSCMainGUI_OutputFcn, ...
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

% --- Executes just before DSCMainGUI is made visible.
function DSCMainGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DSCMainGUI (see VARARGIN)

% Choose default command line output for DSCMainGUI
handles.output = hObject;

% The script of DCEInit is needed for usage of Gilad's and Guy's functions:
% DCEInit;

axes(handles.DSC_curve_fig);
curr_folder=pwd;
imshow([curr_folder,'\DSC_curve.JPG']);

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

% set some strings and values in the GUI:
set(handles.src_folder_str_to_display,'String',[]);
set(handles.processing_status_display,'String',[]);
% Update handles structure
guidata(hObject, handles);



% UIWAIT makes DSCMainGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DSCMainGUI_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on button press in checkbox_PeakSaturation.
function checkbox_PeakSaturation_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PeakSaturation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.filter_PeakSaturation_enable=get(hObject,'Value');
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
display(' ');
display('-------Start pre-processing . . .-------');
DSCCoregFolder=['DSCMainCoreged\'];
dest_dir_str=get(handles.dest_folder_string,'String');
DSCCoregP=[dest_dir_str,filesep,DSCCoregFolder];
% Co-registration
% if '1' - Do the Coregistration if needed or if "Force"==1
if ~(exist(DSCCoregP)&&isdir(DSCCoregP)&&length(dir([DSCCoregP,'\Coreged*.nii']))==length(dir([handles.Nii_src_folder,'\vol*.nii'])) ) || get(handles.checkbox_CoRegistration,'Value')
    Nii_src_folder=handles.Nii_src_folder;
    [DSCCrgDir DSCCrgFNs] = DSCCoreg(DSCCoregP,Nii_src_folder);
    WorkingP=DSCCrgDir;
    display(['Done Co-registration.']);
    display(['Coreged images are in folder: ',WorkingP]);
else
    WorkingP=DSCCoregP;
    display('No Co-registration was needed or forced');
    display(['Images are in folder: ',WorkingP]);
end

% Arrange (coreged) data in 4D, and calc mean
DDSC=dir([WorkingP 'Coreged*vol*.nii']);
DSCFNs=strcat(WorkingP,{DDSC.name})';
data4D_init_temp=loadCoregedNiftis(DSCFNs);
data4D_init=flipdim(permute(data4D_init_temp,[2 1 3 4]),1);
data4DmeanTime=mean(data4D_init,4);
MeanFN=[WorkingP 'DSCMean.nii'];
Raw2Nii(data4DmeanTime,MeanFN,'float32',DSCFNs{1});
handles.data4D_init=data4D_init;
display('Arranging data in 4D...');

% Brain extraction:
TH_BrainExtraction=str2double(get(handles.Brain_TH_textbox,'String'));
[BetStrippedMeanFN, BetMaskFN]=mcbet2(MeanFN,1,TH_BrainExtraction);
ManMaskFN=[WorkingP 'Manual_BrainMask.nii'];
if(exist(ManMaskFN,'file'))
    %AddToLog(WorkingP,'a_2caaaaa','Using Manual_BrainMask');
%     AddToLog(WorkingP,'a_2caaaaa','Using Manual.BrainMask');
    BrainMask=loadniidata(ManMaskFN)>0;
else
    BrainMask=loadniidata(BetMaskFN)>0;
end
data4D_brain_mask=data4D_init.*repmat(BrainMask,[1 1 1 size(data4D_init,4)]);
display(' ');
display('Creating brain mask.');

% Filters for pre-process
% checkbox_LowAtZero_Callback(hObject, eventdata, handles);
% checkbox_LowAverageBaseline_Callback(hObject, eventdata, handles);
% checkbox_NoBolusPeak_Callback(hObject, eventdata, handles);
% checkbox_PeakOutOfBolusRange_Callback(hObject, eventdata, handles);
% checkbox_BigFluctuations_Callback(hObject, eventdata, handles);
handles=guidata(hObject);
handles.BrainMask=BrainMask;
enable2=get(handles.checkbox_LowAtZero,'Value');%    handles.filter_LowAtZero_enable;
enable3=get(handles.checkbox_LowAverageBaseline,'Value');
enable4=get(handles.checkbox_NoBolusPeak,'Value');
enable5=get(handles.checkbox_PeakSaturation,'Value');
enable6=get(handles.checkbox_BigFluctuations,'Value');
enable7=get(handles.checkbox_ZeroValues,'Value');
enable8=get(handles.checkbox_LowSteadyState,'Value');
filter_is_on=[enable2 enable3 enable4 enable5 enable6 enable7 enable8];
display('Chosen filters will be attempted:');


% Apply the filters, not passing the new data4D yet
% data4D_init=handles.data4D_init;
Header=handles.header;
[Mask,data4D_filtered,Filter_per_type,Filter_total] = data4D_mask_filtering(data4D_brain_mask,BrainMask,filter_is_on);
num_voxels_total=numel(Mask);
num_voxels_total_brain=numel(find(BrainMask));
data_voxels_row_inds=Filter_total.data_voxels_row_col_slice(:,1);
data_voxels_col_inds=Filter_total.data_voxels_row_col_slice(:,2);
data_voxels_slice_inds=Filter_total.data_voxels_row_col_slice(:,3);
handles.data4D_filtered=data4D_filtered;
handles.Mask=Mask;
handles.data_voxels_row_col_slice=Filter_total.data_voxels_row_col_slice;
handles.data_voxels_indices=Filter_total.data_voxels_indices;

%Showing filter to user, then passing the data4d to the next step:
setappdata(0,'hMainGUI',gcf);
setappdata(0,'data4D_init',data4D_init);
setappdata(0,'data_voxels_row_col_slice',Filter_total.data_voxels_row_col_slice);
setappdata(0,'Mask',Mask);
setappdata(0,'BrainMask',BrainMask);
% guidata(hObject,handles);
display('Showing total mask to user...');
f=show_filters_GUI;
setappdata(f,'data4D_filtered',data4D_filtered);

% f=gcf;
uiwait(f);
% data_show_filter_GUI = guidata(f);
% apply_filter=data_show_filter_GUI.apply_filter;
apply_filter=getappdata(0,'apply_filter');
close(f);
% if user pressed "back", the filter will not apply, return to main GUI to
% do it again
if ~apply_filter
    set(handles.processing_status_display,'String','Waiting for new user filter-choosing');
    display('Waiting For new user filter-choosing');
    return
end

guidata(hObject,handles);

setappdata(0,'hMainGui',gcf);
display(['Applying filters.' , num2str(double(round(100*length(Filter_total.data_voxels_indices)/num_voxels_total_brain*10)/10)) ,'% of brain voxels passed mask.']);


% Find bolus start, compute baseline  (code taken from Gilad)

disp('Rough estimation of bolus time');
% We mask again for all slices with value bigger than minimum (this time for all time periods and not just the first)
% DCE2D=Reshape4d22d(DCE4D,MskMinSignal);
DSC2D=reshape(data4D_filtered,size(data4D_filtered,1)*size(data4D_filtered,2)*size(data4D_filtered,3),size(data4D_filtered,4));
DSC2D_data=DSC2D(Filter_total.data_voxels_indices,:);

% Get the median of each of the 3d images for the entire time slots 
MedTC=median(DSC2D_data,1);
% [~, optmixture] = GaussianMixture(MedTC', 3, 0,false);
% [Z Grouped]=max(optmixture.pnk,[],2);
% if(max(Grouped)==1)
%     Smoothed=conv2(MedTC,ones(1,ceil(nVols/10)),'same');
%     [a BolusStart]=max(Smoothed);
% else
%     Grouped(end)=Grouped(1)+1;
%     BolusStart=find(Grouped~=Grouped(1),1);end
% Second method
% TwoMinTimePoint=floor(2/TimeBetweenDCEVolsMin);
Ps=zeros(1,numel(MedTC))+2;
% We use the t-test to get the biggest probability that the distribution of the sample
% is diffrent than the rest of the test ( -> smallest Ps value)
for i=3:min(numel(MedTC)-2) %Take the minimum out of 2 minutes frame to 2 frames before the end
    [h Ps(i)]=ttest2(MedTC(1:i),MedTC((i+1):end),[],[],'unequal');
end
mLPs=-log(Ps);
% figure;plot(1:numel(MedTC),MedTC,'b',1:numel(MedTC),mLPs.*(max(MedTC)-min(MedTC))./(max(mLPs)-min(mLPs))+min(MedTC),'r')
[Tmp, BolusStart]=max(mLPs);

handles.BolusStart=BolusStart;
% ASK GILAD - why did he add +1?
% BolusStart=BolusStart+1;
% BolusStartMin=(BolusStart-1)*TimeBetweenDCEVolsMin;
% BolusStart=find(MedTC>MedTC(1)+20,1);
% The base line is the mean of the first images until the bolus
% Baseline=mean(DCE4D(:,:,:,1:(BolusStart-2)),4);
% BaselineFN=[WorkingP 'Baseline.nii'];
% Raw2Nii(Baseline,BaselineFN,'float32', MeanFN);
% figure(78362);subplot(1,2,2);
% plot(MedTC); hold on;plot([BolusStart BolusStart],[min(MedTC) max(MedTC)],'r');
% title('Bolus start approximation');



%Show several DSC data time curves on a separate figure:
display('Showing some random DSC time curves...');
num_DSC_curves_to_show=20;
DSC_curve_indices_to_show=ceil(length(Filter_total.data_voxels_indices)*rand([1 num_DSC_curves_to_show]));
DSC_curves_fig=figure;set(DSC_curves_fig,'units','normalized','outerposition',[0 0 1 1]);
min_val_in_figure=Inf;
max_val_in_figure=-Inf;
for ii=1:num_DSC_curves_to_show
   log_curve_to_show=squeeze(log(data4D_init(data_voxels_row_inds(DSC_curve_indices_to_show(ii)),data_voxels_col_inds(DSC_curve_indices_to_show(ii)),data_voxels_slice_inds(DSC_curve_indices_to_show(ii)),3:end)));
   plot(3:size(data4D_init,4),log_curve_to_show);
   hold all;
   if min(log_curve_to_show)<min_val_in_figure
       min_val_in_figure=min(log_curve_to_show);
   end
   if max(log_curve_to_show)>max_val_in_figure
       max_val_in_figure=max(log_curve_to_show);
   end
end
plot([BolusStart BolusStart],[min_val_in_figure-0.1 max_val_in_figure+0.1],'r','LineWidth',2)
grid on;
grid minor;
axis_vec=axis;
set(gca,'xtick',[0:2:size(data4D_init,4)]);
set(gca,'ytick',[0:axis_vec(4):axis_vec(4)])
title('Ln of DSC intensity curves (randomly chosen)');


%Open the bolus properties GUI and use its output:
setappdata(0,'first_bl_sample',str2double(get(handles.first_bl_sample_value,'String')));
setappdata(0,'last_bl_sample',str2double(get(handles.last_bl_sample_value,'String')));
% setappdata(0,'last_sample',str2double(get(handles.last_sample_value,'String')));
setappdata(0,'last_sample',size(data4D_init,4));
fh=bolus_properties_GUI;
uiwait(fh);
hMainGui=getappdata(0,'hMainGui');
first_bl_sample=getappdata(hMainGui,'first_bl_sample');
last_bl_sample=getappdata(hMainGui,'last_bl_sample');
last_sample=getappdata(hMainGui,'last_sample');

handles.first_bl_sample=first_bl_sample;
handles.last_bl_sample=last_bl_sample;
handles.last_sample=last_sample;

set(handles.first_bl_sample_value,'String',num2str(first_bl_sample));
set(handles.last_bl_sample_value,'String',num2str(last_bl_sample));
set(handles.last_sample_value,'String',num2str(last_sample));

% close(fh);
close(DSC_curves_fig);
display('Got time curve parameters from user.');

% convert the filtered 4D data to concentration c(t):   %%%%%%%%%%%%
% name_str=Header.hdr.hk.db_name;
display('Converting from DSC time curve to concentration c(t)...');
name_str=Header.hdr.hist.descrip;
TE=str2double(name_str(strfind(name_str,'TE')+3:strfind(name_str,'TE')+4));
if isnan(TE)
    TE=30;
elseif TE<0 || TE>1000
    error('TE has an invalid value: ',TE);
end
update_bolus_properties_Callback(hObject, eventdata, handles);
handles=guidata(hObject);
handles.data_average_bl=mean(data4D_filtered(:,:,:,first_bl_sample:last_bl_sample),4);
baseline_edges=[handles.first_bl_sample handles.last_bl_sample];
% [concentration_4D Mask] = Intens2concentration_4D(data4D_filtered,baseline_edges,TE,Mask);
concentration_4D=handles.concentration_4D;
Mask=handles.Mask;
handles.BolusStart=BolusStart;

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
data4D_rgb=repmat(data4D_filtered,[1 1 1 1 3]);
axes(handles.data_window);
set(handles.time_sample_slider,'Enable','on');
set(handles.time_sample_slider,'Max',last_sample);
set(handles.slice_slider,'Enable','on');
set(handles.slice_slider,'Max',size(data4D_filtered,3));
SliderStep=1/(size(data4D_filtered,3)-1);
set(handles.slice_slider,'SliderStep',[SliderStep SliderStep]);
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
handles.data4D_init=data4D_init;

guidata(hObject,handles);

% display "done" message:
set(handles.processing_status_display,'String',['Done. ',num2str(double(round(100*length(Filter_total.data_voxels_indices)/num_voxels_total_brain*10)/10)) ,'% of brain voxels passed mask.']);
display('Done pre-processing!');
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

% --- Executes on button press in Dicom_button
function Dicom_button_Callback(hObject, eventdata, handles)
% hObject    handle to nii_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val=get(hObject,'value');
if val==get(hObject,'Max')
    handles.data_file_type='dicom';
else
    handles.data_file_type='nii';
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
% nii_button_Callback(hObject, eventdata, handles);
% Dicom_button_Callback(hObject, eventdata, handles);
set(handles.get_data_text,'String','');
if get(handles.Dicom_button,'value') && ~get(handles.nii_button,'value')
    file_type='dicom';
elseif ~get(handles.Dicom_button,'value') && get(handles.nii_button,'value')
    file_type='nii';
else
    Display('Problem in choosing Dicom or nii');
end
handles=guidata(hObject);
handles.data_file_type=file_type;
if strcmp(file_type,'dicom')
    file_type_ext='dcm';
else 
    file_type_ext=file_type;
end
file_type_str=strcat('*.',file_type_ext);
init_folder_choose='D:\users\chenh\DSC_project\temp\temp25\Nii_from_Dicoms';
[FileName,PathName] = uigetfile({file_type_str},'Select data (Choose the first file, the whole list will be used)',init_folder_choose);
handles.init_folder_choose=init_folder_choose;
% PathName=uigetdir2('D:\users\chenh\DSC_project\temp\temp12\Nii_from_Dicoms','Select folder to Open');
set(handles.src_folder_str_to_display,'String',PathName);
% data_file_path=strcat(PathName,FileName);
% set(handles.src_folder_str_to_display,'String',strcat('Opened file: ',data_file_path));


% if strcmp(handles.data_file_type,'nii')
%     [Out Header A]=loadniidata(data_file_path);
% else
%     % read Dicom is not support yet
% %     Display('Reading other file than nii is not supported yet');
%     dest_dir_str_user=get(handles.dest_folder_string,'String');
%     if ~isempty(dest_dir_str_user)&& ~strcmp(dest_dir_str_user,'0')
%         dest_dir=dest_dir_str_user;
%     else
%         dest_dir=PathName;
%     end
%     
%     Out=gDicom2Nifti(PathName,TRG_FN)
% end
% arrange the image data in the proper orientation:

guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function src_folder_str_to_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to src_folder_str_to_display (see GCBO)
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
name_str=Header.hdr.hist.descrip;
TE=0.001*str2double(name_str(strfind(name_str,'TE')+3:strfind(name_str,'TE')+4));
if isnan(TE)
    TE=0.030;
elseif TE<0 || TE>1000
    error('TE has an invalid value: ',TE);
end
handles.TE=TE;
data4D_filtered=handles.data4D_filtered;
[concentration_4D Mask] = Intens2concentration_4D(data4D_filtered,baseline_edges,last_sample,TE,Mask,handles);

handles.concentration_4D=concentration_4D;
handles.first_bl_sample=first_bl_sample;
handles.last_bl_sample=last_bl_sample;
handles.last_sample=last_sample;
handles.N_time_points_final=last_sample-first_bl_sample+1;
handles.Mask=Mask;
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
if 0<col && col<=size( handles.data4D_filtered,1) && 0<row && row<=size( handles.data4D_filtered,2)
    set(handles.row_text,'String',num2str(row));
    set(handles.col_text,'String',num2str(col));
    axes(handles.plot_AIFs);
    slice=round(get(handles.slice_slider,'Value'));
    concentration_4D=handles.concentration_4D;
    Mask=handles.Mask;
    AIF=handles.AIF;
    c_t_curve_to_show=squeeze(concentration_4D(row,col,slice,:))*Mask(row,col,slice);
    
%     plot_AIF_
    plot_AIF(AIF,c_t_curve_to_show,handles);
    % choose the max of y axis, for a comfortable view (there are 2
    % different options. One is the global maximum, the other is the max of
    % ~90% from the curves.
%     local_max=max(c_t_curve_to_show);
%     if local_max>handles.c_t_majority_max
%         y_axis_max=handles.c_t_global_max;
%     else
%         y_axis_max=handles.c_t_majority_max;
%     end
%     plot(AIF,'r','LineWidth',1.5);
%     hold on;
%     plot(c_t_curve_to_show);
%     hold off;
%     grid on
%     axis([handles.first_bl_sample handles.last_sample handles.c_t_global_min y_axis_max]);
%     set(gca,'xtick',[1:handles.last_bl_sample-1:handles.last_bl_sample]);
% %     axis([0 length(concentration_4D(row,col,slice,:)) c_t_min-5 c_t_max+5]);
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

%     axes(handles.plot_AIFs);
    slice=round(get(handles.slice_slider,'Value'));
    concentration_4D=handles.concentration_4D;
    Mask=handles.Mask;
    curves4aif_num=curves4aif_num+1;  % update the number of curves participating
    inds4aif_mat(curves4aif_num,:)=[row col slice];
    % compute and plot AIF with the new concentraion curve:
    [AIF,AIF_std]=calc_AIF_from_curves(concentration_4D,inds4aif_mat,'average');
    CTC=AIF;
    plot_AIF(AIF,CTC,handles);

%     % choose the max of y axis, for a comfortable view (there are 2
%     % different options. One is the global maximum, the other is the max of
%     % ~90% from the curves.
%     if max(AIF)>handles.c_t_majority_max
%         y_axis_max=handles.c_t_global_max;
%     else
%         y_axis_max=handles.c_t_majority_max;
%     end
%     plot(AIF,'r','Linewidth',1.5);
%     axis([handles.first_bl_sample handles.last_sample handles.c_t_global_min y_axis_max]);    
    

   % mark the chosen voxel in the brain image (data_window)
    axes(handles.data_window);
    data4D_rgb=handles.data4D_rgb;
    %mark voxel in yellow:
    image_max=repmat((max(max(data4D_rgb(:,:,slice,:,2)))),[1 1 1 1 3]);
    yellow_vec=[1 1 0];
%     yellow_value=image_max*repmat(reshape([1 1 0],1,1,1,1,3);
    yellow_mat=image_max.*repmat(reshape(yellow_vec,1,1,1,1,3),[1,1,1,size(data4D_rgb,4)]);
    data4D_rgb(row,col,slice,:,:)=yellow_mat;
    image_to_show=squeeze(data4D_rgb(:,:,slice,time_sample,:));
    imshow(mat2gray(image_to_show));

    % update curves listbox:
%     String_chen={'Row 5';'Row 5'};
    curves4aif_listbox_string=cell(curves4aif_num,1);
    for ii=1:size(inds4aif_mat,1);
        curves4aif_listbox_string{ii}=['Row ',num2str(inds4aif_mat(ii,1)),'    Col ',num2str(inds4aif_mat(ii,2)),'    Slice ',num2str(inds4aif_mat(ii,3))];
        
    end
    set(handles.curves4aif_listbox,'String',curves4aif_listbox_string);
    
    % update values in handles:
    handles.curves4aif_num=curves4aif_num;
    handles.inds4aif_mat=inds4aif_mat;
    handles.curves4aif_listbox_string=curves4aif_listbox_string;
    handles.AIF=AIF;
    handles.data4D_rgb=data4D_rgb;
    handles.time_sample_to_show=time_sample;
    handles.slice_to_show=slice;
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


% --- Executes on button press in delete_selected_listbox_button.
function delete_selected_listbox_button_Callback(hObject, eventdata, handles)
% hObject    handle to delete_selected_listbox_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
items_to_del=get(handles.curves4aif_listbox,'Value');
inds4aif_mat=handles.inds4aif_mat;
curves4aif_num=handles.curves4aif_num;
curves4aif_listbox_string=handles.curves4aif_listbox_string;
concentration_4D=handles.concentration_4D;
data4D_filtered=handles.data4D_filtered;
time_sample_to_show=handles.time_sample_to_show;
slice_to_show=handles.slice_to_show;
data4D_rgb=handles.data4D_rgb;

%return to the orig values in brain image:
axes(handles.data_window);
for ii=1:length(items_to_del)
    inds4aif_to_del=inds4aif_mat(items_to_del(ii),:);
    data_for_voxel=repmat(data4D_filtered(inds4aif_to_del(1),inds4aif_to_del(2),inds4aif_to_del(3),:),[1 1 1 1 3]);
    data4D_rgb(inds4aif_to_del(1),inds4aif_to_del(2),inds4aif_to_del(3),:,:)=data_for_voxel;
end
image_to_show=squeeze(data4D_rgb(:,:,slice_to_show,time_sample_to_show,:));
imshow(mat2gray(image_to_show));

%delete chosen items from the list:
inds4aif_mat(items_to_del,:)=[];
curves4aif_listbox_string(items_to_del)=[];
curves4aif_num=curves4aif_num-length(items_to_del);

%calc new AIF based on current voxels:
[AIF,AIF_std]=calc_AIF_from_curves(concentration_4D,inds4aif_mat,'average');

%plot the AIF:
plot_AIF(AIF,AIF,handles);

%update the String in the listbox:
set(handles.curves4aif_listbox,'String',curves4aif_listbox_string);
set(handles.curves4aif_listbox,'Value',length(curves4aif_listbox_string));

%save updated data to handle:
handles.items_to_del=items_to_del;
handles.inds4aif_mat=inds4aif_mat;
handles.curves4aif_listbox_string=curves4aif_listbox_string;
handles.curves4aif_num=curves4aif_num;
handles.data4D_rgb=data4D_rgb;
guidata(hObject,handles);



% --- Executes on button press in delete_all_listbox_button.
function delete_all_listbox_button_Callback(hObject, eventdata, handles)
% hObject    handle to delete_all_listbox_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curves4aif_listbox_string=handles.curves4aif_listbox_string;
num_elements_listbox=length(curves4aif_listbox_string);
set(handles.curves4aif_listbox,'Value',[1:num_elements_listbox]);
% use the function for deleting selected, on ALL items in listbox:
delete_selected_listbox_button_Callback(hObject, eventdata, handles)



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
display(' ');
display('-----Start produce maps-----');
display('AIF calculated (user chosed manually)');
AIF=handles.AIF;
concentration4D=handles.concentration_4D;
Mask=handles.Mask;
handles.permeability_correction=get(handles.perm_correction_checkbox,'Value');
deconv_methods.sSVD.en=get(handles.checkbox_sSVD,'Value');
deconv_methods.sSVD.th=str2double(get(handles.checkbox_sSVD_th,'String'));
deconv_methods.cSVD.en=get(handles.checkbox_cSVD,'Value');
deconv_methods.cSVD.th=str2double(get(handles.checkbox_cSVD_th,'String'));
deconv_methods.oSVD.en=get(handles.checkbox_oSVD,'Value');
deconv_methods.oSVD.OI=str2double(get(handles.checkbox_oSVD_OI,'String'));
deconv_methods.tikhonov.en=get(handles.checkbox_tikhonov,'Value');
tic;
[CBF,CBV,MTT,K1,K2,TTP]=calc_maps(concentration4D,AIF,Mask,deconv_methods,handles);
maps_calc_time=toc;
maps_path=save_maps(CBF,CBV,MTT,K1,K2,TTP,Mask,deconv_methods,handles);
maps_folder_message_str=['Maps saved to: ',maps_path];
set(handles.maps_folder_message,'Enable','on');
set(handles.maps_folder_message,'String',maps_folder_message_str);
set(handles.open_maps_folder_button,'Enable','On');
handles.maps_path=maps_path;
display('Done calculating maps!');
display(['(Calculation time: ',num2str(maps_calc_time),' [sec])']);
display(['Maps saved to: ' ,maps_path]);
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
handles=guidata(hObject);

Explore(handles.maps_path);

% --- Executes during object creation, after setting all properties.
function maps_folder_message_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maps_folder_message (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in preview_filtering_button.
function preview_filtering_button_Callback(hObject, eventdata, handles)
% hObject    handle to preview_filtering_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

guidata(hObject,handles);


% --- Executes on selection change in curves4aif_listbox.
function curves4aif_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to curves4aif_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns curves4aif_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from curves4aif_listbox
handles=guidata(hObject);
listbox_all_str= get(hObject,'String');
selected_curve4aif_ind=get(hObject,'Value');
listbox_voxel_str=listbox_all_str{selected_curve4aif_ind};
% find row from 
row_digit1=listbox_voxel_str(strfind(listbox_voxel_str,'Row')+4);
row_digit2=listbox_voxel_str(strfind(listbox_voxel_str,'Row')+5);
if isnan(str2double(row_digit2))
    row=str2double(row_digit1);
else
    row=str2double([row_digit1,row_digit2]);
end

col_digit1=listbox_voxel_str(strfind(listbox_voxel_str,'Col')+4);
col_digit2=listbox_voxel_str(strfind(listbox_voxel_str,'Col')+5);
if isnan(str2double(col_digit2))
    col=str2double(col_digit1);
else
    col=str2double([col_digit1,col_digit2]);
end

slice_digit1_ind=strfind(listbox_voxel_str,'Slice')+6;
slice_digit1=listbox_voxel_str(slice_digit1_ind);
if slice_digit1_ind<length(listbox_voxel_str)
    slice_digit2=listbox_voxel_str(strfind(listbox_voxel_str,'Slice')+7);
    slice=str2double([slice_digit1,slice_digit2]);
else
    slice=str2double(slice_digit1);
end

concentration_4D=handles.concentration_4D;
c_t_curve_to_show=squeeze(concentration_4D(row,col,slice,:));
AIF=handles.AIF;
%  plot_AIF and chosen CTC curve
plot_AIF(AIF,c_t_curve_to_show,handles);
% choose the max of y axis, for a comfortable view (there are 2
% different options. One is the global maximum, the other is the max of
% ~90% from the curves.
% local_max=max(c_t_curve_to_show);
% if local_max>handles.c_t_majority_max
%     y_axis_max=handles.c_t_global_max;
% else
%     y_axis_max=handles.c_t_majority_max;
% end
% plot(AIF,'r','LineWidth',1.5);
% hold on;
% plot(c_t_curve_to_show);
% hold off;
% grid on
% axis([handles.first_bl_sample handles.last_sample handles.c_t_global_min y_axis_max]);
% set(gca,'xtick',[1:handles.last_bl_sample-1:handles.last_bl_sample]);
% %     axis([0 length(concentration_4D(row,col,


% --- Executes on button press in choose_dest_folder_button.
function choose_dest_folder_button_Callback(hObject, eventdata, handles)
% hObject    handle to choose_dest_folder_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
init_folder_choose=handles.init_folder_choose;
dest_dir=uigetdir(init_folder_choose);
set(handles.dest_folder_string,'String',[dest_dir]);
set(handles.get_data_text,'String','');
handles.dest_dir=dest_dir;
guidata(hObject,handles);

% --- Executes on button press in get_data_pushbutton.
function get_data_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to get_data_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
src_folder=get(handles.src_folder_str_to_display,'String');
dest_folder=get(handles.dest_folder_string,'String');
display(['source folder: ',src_folder]);
display(['dest folder: ',dest_folder]);
Nii_src_folder=[dest_folder,filesep,'Nii_from_Dicoms'];
if ~exist(Nii_src_folder)
    mkdir(Nii_src_folder);
end
handles.Nii_src_folder=Nii_src_folder;
if strcmp(handles.data_file_type,'dicom')
    % transform from dicom to nii:
%     set(handles.get_data_text,'String','Opening Dicoms data');
%     guidata(hObject,handles);
    gDicom2Nifti(src_folder,[Nii_src_folder,filesep 'vol.nii']);
else % nii files
    if ~strcmp(src_folder,Nii_src_folder)
        Dir_src=dir([src_folder,'\*.nii']);
        for ii=1:length(Dir_src)
           copyfile([src_folder,filesep,Dir_src(ii).name],[Nii_src_folder,filesep,Dir_src(ii).name]);
        end
    end
end
% Get the header (of one of the nii created before)
% find the exact string of the nifti files, and if begins in 0000 or 0001:
nii_files_list=dir(Nii_src_folder);
nii_files_list_str=[nii_files_list.name];
vol_inds=strfind(nii_files_list_str,'vol');
first_vol_ind=vol_inds(1);
second_vol_ind=vol_inds(2);
first_vol_str=nii_files_list_str(first_vol_ind:second_vol_ind-1);
[Out Header A]=loadniidata([Nii_src_folder,filesep,first_vol_str]);
% voxels_data_4D=flipdim(permute(Out,[2 1 3 4]),1);
% 
% handles.data4D_init=voxels_data_4D;
handles.header=Header;
% handles.data_file_path=data_file_path;
% guidata(hObject, handles);
% % Set the values of the data window sliders according to the data:
% N_slices=size(voxels_data_4D,3);
% N_time_points_init=size(voxels_data_4D,4);
% set(handles.slice_slider,'Enable','on');
% set(handles.slice_slider,'Max',N_slices);
% set(handles.slice_slider,'Min',1);
% set(handles.slice_slider,'SliderStep',[1 1]/(N_slices-1));
% set(handles.slice_slider,'Value',1);
% set(handles.time_sample_slider,'Enable','on');
% set(handles.time_sample_slider,'Max',N_time_points_init);
% set(handles.time_sample_slider,'Min',1);
% set(handles.time_sample_slider,'SliderStep',[1 1]/(N_time_points_init-1));
% set(handles.time_sample_slider,'Value',1);
set(handles.get_data_text,'String','Read data successfully');
display('Read data successfully');
guidata(hObject,handles);


% --- Executes on button press in checkbox_LowSteadyState.
function checkbox_LowSteadyState_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_LowSteadyState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_LowSteadyState


% --- Executes on button press in DCE_init_pushbutton.
function DCE_init_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to DCE_init_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
curr_fold=pwd;
cd ..\..
DSCInit;
set(handles.DSCInit_text,'String','DSCInit finished successfully!'); 
cd(curr_fold);
guidata(hObject,handles);



function Brain_TH_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to Brain_TH_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Brain_TH_textbox as text
%        str2double(get(hObject,'String')) returns contents of Brain_TH_textbox as a double


% --- Executes during object creation, after setting all properties.
function Brain_TH_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Brain_TH_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in perm_correction_checkbox.
function perm_correction_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to perm_correction_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of perm_correction_checkbox


% --- Executes on button press in calcTTP_checkbox.
function calcTTP_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to calcTTP_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calcTTP_checkbox


% --- Executes on button press in AIF_from_file_pushbutton.
function AIF_from_file_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to AIF_from_file_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guidata(hObject);
file_type_str={'*.txt','*.mat'};
init_folder_choose=handles.init_folder_choose;
[FileName,PathName] = uigetfile(file_type_str,'Choose AIF curve file',init_folder_choose);
AIFcurve_text_array=textread([PathName FileName],'%s');
[Lia,Loc]=ismember('bolusinfo',AIFcurve_text_array);
if Lia==1
    handles.first_bl_sample=str2double(AIFcurve_text_array(Loc+2));
    set(handles.first_bl_sample_value,'String',AIFcurve_text_array(Loc+2));
    handles.last_bl_sample=str2double(AIFcurve_text_array(Loc+3));
    set(handles.last_bl_sample_value,'String',AIFcurve_text_array(Loc+3));
    handles.last_sample=str2double(AIFcurve_text_array(Loc+4));
    set(handles.last_sample_value,'String',AIFcurve_text_array(Loc+4));
    update_bolus_properties_Callback(hObject, eventdata, handles);
    handles=guidata(hObject);
    AIF=str2double(AIFcurve_text_array(Loc+5:end));
    %in the file, the AIF is 1000 times smaller and need to be normalized:
    AIF=AIF*1000;
    handles.AIF=AIF;
    plot_AIF(AIF,AIF,handles);
else
   display('Cannot reaf AIF from file. The text file is not in the right format');
   set(handles.AIF_from_file_textbox,'String','Error reading AIF from file');
end
guidata(hObject,handles);
