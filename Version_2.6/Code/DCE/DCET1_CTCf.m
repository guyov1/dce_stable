% The SPGR equation:
% S(t) = M0 * sin(FA)*( 1-E1(t) ) / ( 1-E1(t)cos(FA) )
% where E1=exp(-TR/T1(t))
% So, given M0 from relaxometry (and TR, FA, S), we can get T1(t)
% The contrast agent concentration equation
% 1/T1(t) = 1/(T1base)+C(t)
% So given T1base, we can get from T1(t) the Concentration Time Curve C(t)

function DCET1_CTCf(DCEInfo, Additional_T1_Maps_Time_Stamps, WorkingP,DoN3,DoGlobal,DoDeviations,CalcForce,Options)

% If using super-resolution, get the relevant parameters value
UnderSampling=Options.SubSampling;
if(UnderSampling==1)
    USStr='';
else
    USStr=['_' num2str(UnderSampling)];
end

% .mat file name for that stage
CTCFN=[WorkingP 'AfterCTC' USStr '.mat'];

% If the .mat file of this stage alreay exists, return.
if(exist(CTCFN,'file') && ~CalcForce)
    %     TmpWorkingP=WorkingP;
    %     load(CTCFN);
    %     WorkingP=TmpWorkingP;
    disp('DCE_CTC Already computed');
    return;
end

AddToLog(WorkingP,'c_00','\\subsection*{Extracting the CTCs}');

% Starting CTC calculation message
disp('DCE_CTC Starting');

% Get the mean DCE map
MeanFN=[WorkingP 'DCEMean.nii'];
%% Load all previous stages
PrepareFN=[WorkingP 'AfterPrepare4D.mat'];
load(PrepareFN,'Baseline','Msk','DCE4D','TimeBetweenDCEVols','BolusStart','BrainMask');
if(isfield(Options,'TimeMultiplier'))
    TimeBetweenDCEVolsFinal=TimeBetweenDCEVols*Options.TimeMultiplier;
end
if(isfield(Options,'TimeMultiplier') && Options.TimeMultiplier<0)
    TimeBetweenDCEVolsFinal=-Options.TimeMultiplier;
end
%% Update 4D map according to additional DESPOT1 data we have

% Get all T1 maps we have
All_T1_Maps = dir( [WorkingP 'Relaxometry' filesep 'Set*']);

% Get number of T1 maps
Num_Of_T1_Maps = numel(All_T1_Maps);

% If we have more than 1, use it and create time vector
WhichSetForBase=1;
if (Num_Of_T1_Maps>1)
    
    % Initiate time vector
    T1_Maps_Time_Vector = zeros(1,Num_Of_T1_Maps-1);
    % Initiate date vector
    T1_Maps_Time_Date      = zeros(1,Num_Of_T1_Maps-1);
    % Initiate T1 additional maps
    Additional_T1_Maps = NaN([size(Baseline) Num_Of_T1_Maps]);
    CAdditional_T1_Maps = Additional_T1_Maps;
    OAdditional_T1_Maps = Additional_T1_Maps;
    
    %%
    SingleAngleIdxs=[];
    % Ignore the first map which is used for T1-base calculation
    for map_idx = 1:Num_Of_T1_Maps
        dir_to_look = All_T1_Maps(map_idx);
        dir_to_look.name=['Set_Num_' num2str(map_idx)];
        dir_to_look_name = [WorkingP 'Relaxometry' filesep dir_to_look.name filesep 'Middle_Time_Stamp.mat'];
        %         GGG Correct this - single angle
        if(~exist(dir_to_look_name,'file'))
            continue;
        end
        Mid_Time_Stamp = load(dir_to_look_name) ;
        Mid_Time_Stamp = struct2array(Mid_Time_Stamp);
        Mid_Time_Stamp = Mid_Time_Stamp{1};
        % Seperate date and time
        ind_date_time = findstr(Mid_Time_Stamp,'_');
        Mid_Time_Stamp_Date = Mid_Time_Stamp(1:ind_date_time-1);
        Mid_Time_Stamp_Time = Mid_Time_Stamp(ind_date_time+1:end);
        
        % Convert time to seconds
        hours       = str2double(Mid_Time_Stamp_Time(1:2));
        minutes  = str2double(Mid_Time_Stamp_Time(3:4));
        seconds  = str2double(Mid_Time_Stamp_Time(5:6));
        time_seconds = hours*3600 + minutes*60 + seconds;
        
        % Update time and date vectors
        T1_Maps_Time_Vector(map_idx) = time_seconds;
        T1_Maps_Time_Date(map_idx)      = str2double(Mid_Time_Stamp_Date);
        
        % Update T1 nifti maps
        T1_Path = [WorkingP 'Relaxometry' filesep dir_to_look.name filesep];
        
        T1_Map_Path=[T1_Path 'Seg_ForSeg'  filesep 'ForSeg.nii'];
        CT1_Map_Path=[T1_Path 'Seg_ForSeg'  filesep 'mForSeg.nii'];
        OT1_Map_Path=[T1_Path 'T13DOFA.nii'];
        
        %        if (DoN3+DoGlobal+1 == 3)
        %            if (DoDeviations+1 == 2)
        %
        %                T1_Map_Path=[T1_Path 'T13DNFA_N3k.nii'];
        %            else
        %                T1_Map_Path=[T1_Path 'T13DOFA_N3k.nii'];
        %            end
        %        else
        %            if (DoDeviations+1 == 2)
        %                T1_Map_Path=[T1_Path 'T13DNFA.nii'];
        %
        %            else
        %                T1_Map_Path=[T1_Path 'T13DOFA.nii'];
        %            end
        %
        %        end
        DOrig=dir([T1_Path 'OrigFA_*.nii']);
        % g corrected no segmentation on T1 now - if(exist(T1_Map_Path,'file'))
        if(numel(DOrig)>1)
%             Additional_T1_Maps(:,:,:,map_idx) = loadniidata(T1_Map_Path);
%             CAdditional_T1_Maps(:,:,:,map_idx) = loadniidata(CT1_Map_Path);
%             OAdditional_T1_Maps(:,:,:,map_idx) = loadniidata(OT1_Map_Path);
        else % Single Angle
            SingleAngleIdxs=[SingleAngleIdxs map_idx];
        end
    end
    B1Maps=sqrt(CAdditional_T1_Maps./Additional_T1_Maps);
    %%
    % Get the main data time and date
    num_volumes = size(DCE4D,4);
    num_slices = size(DCE4D,3);
    Main_First_Time_Stamp = DCEInfo.SeriesDateTime;
    % Seperate date and time
    ind_date_time = findstr(Main_First_Time_Stamp,'_');
    Main_First_Time_Stamp_Date = Main_First_Time_Stamp(1:ind_date_time-1);
    Main_First_Time_Stamp_Date = str2double(Main_First_Time_Stamp_Date);
    Main_First_Time_Stamp_Time = Main_First_Time_Stamp(ind_date_time+1:end);
    % Convert time to seconds
    hours       = str2double(Main_First_Time_Stamp_Time(1:2));
    minutes  = str2double(Main_First_Time_Stamp_Time(3:4));
    seconds  = str2double(Main_First_Time_Stamp_Time(5:6));
    time_seconds_main = hours*3600 + minutes*60 + seconds;
    
    
    % Change the time vector to be relative to the first time index ( +1 so index won't be zero)
    time_seconds_main = time_seconds_main - T1_Maps_Time_Vector(1)  + 1;
    
    T1_Maps_Time_Vector = T1_Maps_Time_Vector - T1_Maps_Time_Vector(1) + 1;
    time_seconds_main = time_seconds_main - T1_Maps_Time_Vector(1)  + 1;
    % If the date is different, the test crossed midnight -> seconds calculation is not correct
    GoodAnglesF=find(T1_Maps_Time_Date>0);
    if (length(unique([T1_Maps_Time_Date(GoodAnglesF) Main_First_Time_Stamp_Date])) ~= 1)
        error(' Error! - The dates of the T1 maps are different! Probably did it around midnight!');
    end
    
    
    % Create the time vector for the main acquisition
    % I assumed the acquisitions are uniformly distributed from the first
    % main image time to the first image after the main time
    all_bigger_than_main_time = T1_Maps_Time_Vector(T1_Maps_Time_Vector>time_seconds_main);
    
    % Check if there are no T1 maps after main
    no_maps_after_main = isempty(all_bigger_than_main_time);
    if (no_maps_after_main)
        
        % If there are no maps after, create a time vector from first main
        % time to last main time with TimeBetweenDCEVols seconds apart.
        num_volumes = size(DCE4D,4);
        last_volume_time = time_seconds_main + (num_volumes-1)*TimeBetweenDCEVolsFinal;
        time_diff = TimeBetweenDCEVolsFinal;
        Main_Vector_Time = round(time_seconds_main:time_diff:last_volume_time);
        
    else % There are maps after main
        
        first_after_main_time = all_bigger_than_main_time(1);
        Main_Vector_Time = round(linspace(time_seconds_main,first_after_main_time,num_volumes+1));
        Main_Vector_Time = Main_Vector_Time(1:end-1);
        
    end
    
    % Create the complete time vector (used in FindPKBATgAIFMurase.m)
    Complete_Time_Vector   = sort([T1_Maps_Time_Vector Main_Vector_Time]);
    
    Additional_T1_Maps_Time_Diff=datenum(Additional_T1_Maps_Time_Stamps, 'yyyymmdd_HHMMSS')-datenum(DCEInfo.SeriesDateTime,'yyyymmdd_HHMMSS')';
    a=load([WorkingP 'Relax.mat']);
    num_before=sum(Additional_T1_Maps_Time_Diff(a.index_new_sets)<0);
    Additional_T1_Maps_Time_Diff_Sets=Additional_T1_Maps_Time_Diff(a.index_new_sets);
    
    WhichSetForBase=max(setdiff(find(Additional_T1_Maps_Time_Diff_Sets<0),SingleAngleIdxs));
end
WhichSetForBase=1;
%%
% Get the T1s data from the first set calculated
RWorkingP=[WorkingP 'Relaxometry' filesep 'Set_Num_' num2str(WhichSetForBase) filesep];
FMaskFN=[RWorkingP 'FBrainMsk.nii'];

% Get the T1 base map (the second one is FA corrected - if exists)
T1Res{1,1,1}=[RWorkingP 'T13DOFA.nii'];
T1Res{1,2,1}=[RWorkingP 'T13DNFA.nii'];

% ASK GILAD - Didn't he run the b1 cleaning in DCET1_RelaxForSubjectf?
% ASNWER  - He used it just to get the file name (T1Res{1,1,2}).
%                      If it was already calculated, it won't be calculated again(using "false" option).
% if(DoN3)
%     T1Res{1,1,2}=Iterative_n3_b1_clean_mask_FN(3,T1Res{1,1,1},FMaskFN,false,[]);
%     if(DoDeviations)
%         T1Res{1,2,2}=Iterative_n3_b1_clean_mask_FN(3,T1Res{1,2,1},FMaskFN,false,[]);
%     end
% end
%
% % Get the N3 corrected T1 maps
% T1Res{1,1,3}=[RWorkingP 'T13DOFA_N3k.nii'];
% T1Res{1,2,3}=[RWorkingP 'T13DNFA_N3k.nii'];

T1CleanedFN=[RWorkingP 'Seg_ForSeg'  filesep 'mForSeg.nii'];
% T1UnCleanedFN=[RWorkingP 'Seg_ForSeg'  filesep 'ForSeg.nii'];
T1UnCleanedFN=[RWorkingP 'T13DNFA.nii'];
% RelaxFN=T1CleanedFN;
RelaxFN=T1UnCleanedFN;
% Load .nifti data of the T1 map
CT1=loadniidata(RelaxFN);

% Find B1 by WM ref
D = dir([WorkingP 'RefT1_*.nii']);
if(isempty(D))
    D=dir([WorkingP 'RefAuto_Base_WM*.nii']);
end

RefFN=[WorkingP D(1).name];
AddToLog(WorkingP,'c_20Norm',['T1 Normalization ref file: ' strrep(D(1).name,'_','-')]);
RefTrgVal=str2num(RefFN(find(RefFN=='_',1,'last')+1:end-4));
RefVol=loadniidata(RefFN)>0;
D2 = dir([WorkingP 'WMExMask.nii']);
if(~isempty(D2))
    AddToLog(WorkingP,'c_20NormEx','Using WM exclusion mask!');
    ExMask=loadniidata([WorkingP 'WMExMask.nii'])>0;
    RefVol=RefVol & ~ExMask;
end
RefVolX=RefVol & CT1>100 & Baseline>100;
load(PrepareFN,'BadSlicesF2');
% If removing the bad slices, leave us with no remaining voxels, keep them
TmpRefVol = RefVolX;
TmpRefVol(:,:,BadSlicesF2) = false;
if sum(sum(sum(TmpRefVol))) ~= 0
    RefVolX(:,:,BadSlicesF2) = false;
end  
clear TmpRefVol;

%%
FMaskFN=[RWorkingP 'FBrainMsk.nii'];
% Load the brain mask
FBrainMask=loadniidata(FMaskFN)>0;
GoodSlices=setdiff(1:size(FBrainMask,3),BadSlicesF2);
% ASK GILAD - what is the purpose of the following?
if(Options.UseN3OnT1)
    %     PseudoB1=sqrt(double(loadniidata(T1Res{1,DoDeviations+1,2}))./double(loadniidata(T1Res{1,DoDeviations+1,1})));
    %     UncleanedT1=loadniidata(T1Res{1,1+Options.ExtractFAs,1});
    UncleanedT1=loadniidata(T1UnCleanedFN);
%     PseudoB1=sqrt(CT1./UncleanedT1);
    ShowSurfFit=false;
    NonLinTransformPower=1;
    NonLinTransform=@(x) log(x.^(1/NonLinTransformPower));
    rNonLinTransform=@(x) exp(x).^NonLinTransformPower;

    FindB1gWM_LS;
    %%
%     if(ShowSurfFit)
    gfig(5001);clf;
    for i=1:numel(GoodSlices)
            CurSli=GoodSlices(i);
%         CurSli=11;
        F=find(I{3}==CurSli);
        FF=find(IF{3}==CurSli);
        [X Y]=meshgrid(1:size(RefVol,1),1:size(RefVol,2));
        X(:)=IF{1}(FF);
        Y(:)=IF{2}(FF);
        gsubplot(numel(GoodSlices),i);
        CurFld=squeeze(NewField(:,:,CurSli));
        surf(X,Y,CurFld,'EdgeColor','none');
        colormap hsv
        alpha(.7)
        hold on;
        if(isempty(F))
%             plot3(I{1}(F), I{2}(F),rNonLinTransform(SolVec(F))/RefTrgVal,'ko');
%             ax=axis;
%             axis([min(I{1}(F))-10 max(I{1}(F))+10 min(I{2}(F))-10 max(I{1}(F))+10 ax(5:6)]);
        else
            plot3(I{1}(F), I{2}(F),rNonLinTransform(SolVec(F))/RefTrgVal,'ko');
            ax=axis;
            axis([min(I{1}(F))-10 max(I{1}(F))+10 min(I{2}(F))-10 max(I{2}(F))+10 ax(5:6)]);
        end
        title(CurSli);
        zlabel([mean(CurFld(:)) std(CurFld(:))]);
    end
    %%
    MaximizeSaveCloseAndAddToLog(WorkingP,'B1gWM','yc_30B1gWMrRMS','$B_1$ given WM.');
    %%
    B1FN=[RWorkingP 'N3B1.nii'];
    Raw2Nii(PseudoB1,B1FN,'float32', MeanFN);
end
%% Setting threshold values

% ASK GILAD - how did he determine those values?
% ASNWER - He uses those values to see if the T1 calculation failed.
%                    He checks whether 1000 voxels or more has a value greater than T1 = 4000.
%                    He used T1=1200 before for the iterations for calculating the REAL T1_base to converge.
%                    T1=40,000 before was just a really high limit (so we won't realty limit the possible results).

% Maximum value of T1 in image
T1Thresh=10000;

% Maximum number of values exceeding the T1 threshold
NaTThresh=1000;

% Go over each slice
for i=1:size(CT1,3)
    CurMsk=squeeze(Msk(:,:,i));
    CurCT1=squeeze(CT1(:,:,i));
    nAT(i)=sumn(CurCT1(CurMsk)>T1Thresh);
end
BadSlicesAgain=find(nAT>NaTThresh);
% If the number of pixels in slice, exceeding T1Thresh is bigger than NaTThresh, mark the slice as NaN.
CT1(:,:,BadSlicesAgain)=NaN;
%% Calculate M0

load(PrepareFN,'nVols','Min2Timestep','BrainMask');

% Get FA and TR
DCETR=DCEInfo.RepetitionTime;
DCEFA=DCEInfo.FlipAngle;

% Data saved at the end of DCET1_RelaxForSubjectf.m
MatFN=[RWorkingP 'NFARes.mat'];
load(MatFN);
%%
% ASK GILAD - I don't understand the B1 calculations
% ASNWER       - SHOULD READ THE ARTICLE AND CHECK THAT.
if(Options.UseN3OnT1)
    B1=loadniidata(B1FN);
else
    B1=CT1*0+1;
end

S=sort(CT1(RefVolX));
WhichPercent=0.5;
RefCurVal=S(max(1,floor(numel(S)*WhichPercent)));
% RefCurVal=median(CT1(RefVolX));
kCoeff=RefTrgVal/RefCurVal;
AddToLog(WorkingP,'c_1kc',['kCoeff:' num2str(kCoeff)]);
FinalT1=CT1*kCoeff;
FinalT1FromDESPOT1=FinalT1;

M0FN=[RWorkingP 'PD3DNFA.nii'];
M03DNoB1Nok=loadniidata(M0FN);
if(~Options.UseN3OnT1)
    kCoeff=1;
end
% kCoeff=1;
B1Nok=B1;
B1=B1*sqrt(kCoeff);
FinalT1FN=[WorkingP 'FinalT1.nii'];
Raw2Nii(FinalT1,FinalT1FN,'float32',MeanFN);

WMD=dir([WorkingP 'RefT1_*WM*.nii']);
if(isempty(WMD))
    WMD=dir([WorkingP 'Ref*.nii']);
    WMD=WMD(strhas({WMD.name},'WM'));
end
RefWMFN=[WorkingP WMD(1).name];
RefWMVol=loadniidata(RefWMFN)>0;
RefWMVal=median(FinalT1(RefWMVol));
AddToLog(WorkingP,'c_2wm',['RefWMVal:' num2str(RefWMVal)],[],2);

if(Options.ExtractFAs && Options.IncludeMainInT1)
    Tmp =load([WorkingP 'Relaxometry' filesep 'Set_Num_1' filesep 'NFARes.mat']);
    DCEFA=Tmp.ONFAs(end);
    %     DCEFA=ONFAs(end);
end
FA3DNoB1Nok=CT1*0+DCEFA;
FA3D=(FA3DNoB1Nok./B1Nok)/sqrt(kCoeff);
M03D=(M03DNoB1Nok.*B1Nok)*sqrt(kCoeff);
% CT1 is UncleanedT1 * B1^2
%%
% figure(191919);clf;
% subplot(2,2,1);
% imagesc(mritransform(M03DNoB1Nok(GoodRows,GoodCols,MidSli)),[0 25000]);
% subplot(2,2,2);
% imagesc(mritransform(B1(GoodRows,GoodCols,MidSli)),[0.9 1.1]);
% subplot(2,2,3);Tmp=M03DNoB1Nok.*B1;
% imagesc(mritransform(Tmp(GoodRows,GoodCols,MidSli)),[0 22000]);
%%
% M0 =  ( S(t0) * (1-E1(t)cos(FA)) ) / sin(FA)*(1-E1(t))

% The base line is the mean of the first images until the bolus (DCET1_Prepare4Df.m)
DCEToM0=@(T1x,FAx) (Baseline.*(1-exp(-DCETR./T1x).*cosd(FAx)))./((1-exp(-DCETR./T1x)).*sind(FAx));
% DCEM0NoB1Nok=DCEToM0(UncleanedT1,FA3DNoB1Nok);
DCEM0=DCEToM0(FinalT1,FA3D);

% Check if user wanted a single M0
% If so, take the M0 calculated from the first DESPOT
if (Options.Use_Single_M0)
    DCEM0 = loadniidata([WorkingP 'Relaxometry' filesep 'Set_Num_1' filesep 'PD3DNFA.nii']);
end

% DCEM0=(Baseline.*(1-E1.*CosFA))./((1-E1).*SinFA);
DCEM0FN=[WorkingP 'DCE_M0.nii'];
Raw2Nii(DCEM0,DCEM0FN,'float32',MeanFN);

% Get the T1s data from the first set calculated
% DCERelaxP=[WorkingP 'Relaxometry' filesep 'Set_Num_1' filesep];

% ASK GILAD - what is the purpose of the following 2 lines?
% ANSWER       - He wanted to take a little bigger mask than the calculated one.
se=strel('disk',10,8);
DBrainMask=imdilate(FBrainMask,se);

se=strel('disk',4,8);
EBrainMask=imerode(FBrainMask,se);
%%
[Tmp, FAPlus, FAMinus, GoodDiscriminantaPlus, GoodDiscriminantaMinus]=CalcFAgSigT1M0(Baseline(EBrainMask),FinalT1(EBrainMask),M03D(EBrainMask),DCETR,100,1,1);
Tmp2=NaN*Baseline(EBrainMask);
Tmp2(GoodDiscriminantaPlus)=FAPlus;
Tmp2(Tmp2<0)=NaN;
Tmp=Baseline*NaN;
Tmp(EBrainMask)=Tmp2./DCEFA;
Tmp(FinalT1>2000)=NaN;
Tmp(FinalT1<500)=NaN;
Tmp(:,:,BadSlicesF2)=NaN;
BaselineB1Fld=FitNonLinPolyField(2,NonLinTransform,rNonLinTransform,Tmp,isfinite(Tmp));
BaselineFA3D=DCEFA.*BaselineB1Fld;

[Tmp, FAPlus, FAMinus, GoodDiscriminantaPlus, GoodDiscriminantaMinus]=CalcFAgSigT1M0(Baseline(EBrainMask),FinalT1(EBrainMask),DCEM0(EBrainMask),DCETR,100,1,1);
Tmp2=NaN*Baseline(EBrainMask);
Tmp2(GoodDiscriminantaPlus)=FAPlus;
Tmp2(Tmp2<0)=NaN;
Tmp=Baseline*NaN;
Tmp(EBrainMask)=Tmp2./DCEFA;
Tmp(FinalT1>2000)=NaN;
Tmp(FinalT1<500)=NaN;
Tmp(:,:,BadSlicesF2)=NaN;
BaselineB1Fldx=FitNonLinPolyField(2,NonLinTransform,rNonLinTransform,Tmp,isfinite(Tmp));
BaselineFA3Dx=DCEFA.*BaselineB1Fldx;
%%
MidSli=floor(size(DCE4D,3)/2);
RelaxFNx=[WorkingP 'Relax.mat'];
% load(RelaxFNx,'GoodSlices','GoodRows','GoodCols');

Tmp=max(FBrainMask,[],3);
F=find(max(Tmp,[],2));
GoodRows=F(1):F(end);
F=find(max(Tmp,[],1));
GoodCols=F(1):F(end);
%%
figure(99);clf;
[Mn Mx]=FindDR(DCEM0(FBrainMask));
gsubplot(1,3,1,1);
imagesc(mritransform(DCEM0(GoodRows,GoodCols,MidSli)),[0 Mx]);
title('M0 from baseline');
set(gca,'XTick',[],'YTick',[]);
gsubplot(1,3,1,2);
imagesc(mritransform(M03D(GoodRows,GoodCols,MidSli)),[0 Mx]);
title('M0 from DESPOT1');
set(gca,'XTick',[],'YTick',[]);
gsubplot(1,3,1,3);
imagesc(mritransform(M03D(GoodRows,GoodCols,MidSli)./DCEM0(GoodRows,GoodCols,MidSli)),[0.8 1.2])
set(gca,'XTick',[],'YTick',[]);
title('Ratio [0.8 1.2]');
%%
saveas(99,[WorkingP 'M0Compare.png']);
saveas(99,[WorkingP 'M0Compare.fig']);
close(99);
%%
% DCEToT1=@(FAx,M0x) 1./(-log(((Baseline./M0x)-sind(FAx))./((Baseline./M0x).*cosd(FAx)-sind(FAx)))/DCETR);

DCEVolToT1=@(FAx,M0x,Vol) 1./(-log(((Vol./M0x)-sind(FAx))./((Vol./M0x).*cosd(FAx)-sind(FAx)))/DCETR);
% DCEToT1=@(FAx,M0x) DCEVolToT1(FAx,M0x,Baseline);
DCEVolToT1x=@(Vol) DCEVolToT1(FA3D,DCEM0,Vol);
% g corrected for changing TR (additional reps)
DCEVolToT1TR=@(FAx,M0x,Vol,TR) 1./(-log(((Vol./M0x)-sind(FAx))./((Vol./M0x).*cosd(FAx)-sind(FAx)))/TR);
% DCEToT1=@(FAx,M0x,TR) DCEVolToT1TR(FAx,M0x,Baseline,TR);
DCEVolToT1xTR=@(Vol,TR) DCEVolToT1TR(FA3D,DCEM0,Vol,TR);
% g end correction

% T1FromBaseline=DCEToT1((FA3DNoB1Nok./B1Nok)/sqrt(kCoeff),(M03DNoB1Nok.*B1Nok)*sqrt(kCoeff));
% using the M0 from computed from the baseline
% T1FromBaseline=DCEToT1(FA3D,DCEM0);

% using the DESPOT1 M0
% T1FromBaseline=DCEToT1(FA3D,M03D);
% T1FromBaseline(imag(T1FromBaseline)~=0)=0;

% Using the B1 field computed from the baseline
% T1FromBaseline=DCEToT1(BaselineFA3D,M03D);
% T1FromBaseline=DCEToT1(BaselineFA3Dx,M03D);
% T1FromBaseline=DCEToT1(FA3D,DCEM0);
T1FromBaseline=DCEVolToT1x(Baseline);
T1FromBaseline(imag(T1FromBaseline)~=0)=NaN;
FinalT1=T1FromBaseline;
%%
figure(9384);clf;
gsubplot(4,1);
imagesc(mritransform(FinalT1FromDESPOT1(GoodRows,GoodCols,MidSli)),[0 3000]);
set(gca,'XTick',[],'YTick',[]);
colorbar;
title('Final T_1');
gsubplot(4,2);
imagesc(mritransform(T1FromBaseline(GoodRows,GoodCols,MidSli)),[0 3000]);
set(gca,'XTick',[],'YTick',[]);
title('T_1 from baseline');
gsubplot(4,3);
Tmp=FinalT1FromDESPOT1-T1FromBaseline;
imagesc(mritransform(Tmp(GoodRows,GoodCols,MidSli)),[-100 100]);
set(gca,'XTick',[],'YTick',[]);
title('T_1 difference [-100 100]');
gsubplot(4,4);
Xs=-500:500;
Mskx=FinalT1FromDESPOT1<2000 & FinalT1FromDESPOT1>400 & EBrainMask;
Mskx(:,:,BadSlicesF2)=false;
Tmp2=Tmp(Mskx);
bar(Xs,histc(Tmp2,Xs));
title(['Mode: ' num2str(mode(round(Tmp2))) ' Med abs: ' num2str(median(abs(round(Tmp2)))) ]);
%%
saveas(9384,[WorkingP 'T1DCECompare.png']);
saveas(9384,[WorkingP 'T1DCECompare.fig']);
close(9384);
AddToLog(WorkingP,'yc_5',['T1 from DESPOT1 and from baseline of DCE.'],'T1DCECompare.png');
%%
% ShowT1InFirstAndLast=false;
% if(ShowT1InFirstAndLast)
%     CTC2DX=DCEtoCTCf(DBrainMask,DCE4D(:,:,:,1),CT1,FA3D,DCETR, Baseline);
%     T1FromFirstVol=DBrainMask*0;
%     T1FromFirstVol(DBrainMask)=CTC2DX;
%     T1FromFirstVol=1./(T1FromFirstVol+1./CT1);
%     T1FromFirstVol(imag(T1FromFirstVol)~=0)=NaN;
%     CTC2DX=DCEtoCTCf(DBrainMask,DCE4D(:,:,:,end),CT1,FA3D,DCETR, Baseline);
%     T1FromLastVol=DBrainMask*0;
%     T1FromLastVol(DBrainMask)=CTC2DX+1./CT1(DBrainMask);
%     T1FromLastVol=1./T1FromLastVol;
%     T1FromLastVol(imag(T1FromLastVol)~=0)=NaN;
%
%     figure(1002);clf;
%     gsubplot(2,3,1,1);
%     imagesc(mritransform(CT1(:,:,MidSli)),[0 4000]);
%     gsubplot(2,3,1,2);
%     imagesc(mritransform(T1FromFirstVol(:,:,MidSli)),[0 4000]);
%     gsubplot(2,3,1,3);
%     imagesc(mritransform(T1FromLastVol(:,:,MidSli)),[0 4000]);
%     gsubplot(2,3,2,2);
%     Tmp=CT1-T1FromFirstVol;
%     Tmp(~isfinite(Tmp))=0;
%     imagesc(mritransform(Tmp(:,:,MidSli)),[-100 100]);
%     title(median(abs(Tmp(isfinite(Tmp) & Tmp~=0))));
%     gsubplot(2,3,2,3);
%     Tmp=CT1-T1FromLastVol;
%     Tmp(~isfinite(Tmp))=0;
%     imagesc(mritransform(Tmp(:,:,MidSli)),[-100 100]);
% end
%%
% Compute two dimenstional C(t)
NumVols=size(DCE4D,4)-Options.nVolsToRemoveFromEnd;

% New calculation volume by volume
% T1FromBaseline=DCEToT1(BaselineFA3D,M03D);
T12D=NaN([sum(DBrainMask(:)) NumVols]);
for i=1:NumVols
    disp(i);
    CurVol=squeeze(DCE4D(:,:,:,i));
%     CurT1=DCEVolToT1(FA3D,DCEM0,CurVol);
%     CurT1=DCEVolToT1(BaselineFA3D,DCEM0,CurVol);
%     CurT1=DCEVolToT1(BaselineFA3Dx,DCEM0,CurVol);
    CurT1=DCEVolToT1x(CurVol);
    CurT1(imag(CurT1)~=0)=NaN;
    T12D(:,i)=CurT1(DBrainMask);
end

CTC2DBig=(1./T12D)-repmat(1./FinalT1(DBrainMask),[1 NumVols]);
% Check new calculation
% Tmp=abs((CTC2DBig-CTC2DBigx)./CTC2DBig);
% [mean(Tmp(isfinite(Tmp))) median(Tmp(isfinite(Tmp)))]

% Old calculation
% Note - The following is a heavy computation and can freeze the computer for a minute
% CTC2DBig=DCEtoCTCf(DBrainMask,DCE4D(:,:,:,1:(end-Options.nVolsToRemoveFromEnd)),FinalT1,FA3D,DCETR, Baseline);

% TmpAA=CTC2DBig(:,end);
% TmpX=DCEVolToT1(FA3D,DCEM0,squeeze(DCE4D(:,:,:,end)));
% TmpBB=1./TmpX(DBrainMask)-1./FinalT1(DBrainMask);
% OKI=find(isfinite(TmpAA));
% [TmpAA(OKI(1:10)) TmpBB(OKI(1:10))]'
% max([TmpAA(OKI)-TmpBB(OKI)]')
%%
load(PrepareFN,'BadTimePoints');
BadTimePoints=BadTimePoints(BadTimePoints<size(CTC2DBig,2));
for ii=1:numel(BadTimePoints)
    CurPoint=BadTimePoints(ii);
    if(all(~ismember([CurPoint-2 CurPoint-1 CurPoint+1 CurPoint+2],BadTimePoints)) && ( CurPoint+2 <= size(CTC2DBig,2) ))
        CTC2DBig(:,CurPoint)=mean(CTC2DBig(:,[CurPoint-2 CurPoint-1 CurPoint+1 CurPoint+2]),2);
    else
        if(all(~ismember([CurPoint-1 CurPoint+1],BadTimePoints)) && ( CurPoint+1 <= size(CTC2DBig,2) ))
            CTC2DBig(:,CurPoint)=mean(CTC2DBig(:,[CurPoint-1 CurPoint+1]),2);
        else
            if( (all(~ismember([CurPoint-2 CurPoint+2],BadTimePoints)) ) && ( CurPoint+2 <= size(CTC2DBig,2) ))
                CTC2DBig(:,CurPoint)=mean(CTC2DBig(:,[CurPoint-2 CurPoint+2]),2);
            else
                if( (all(~ismember([CurPoint-2 CurPoint-1 CurPoint+2 CurPoint+3],BadTimePoints)) ) && ( CurPoint+3 <= size(CTC2DBig,2) ))
                    CTC2DBig(:,CurPoint)=mean(CTC2DBig(:,[CurPoint-2 CurPoint-1 CurPoint+2 CurPoint+3]),2);
                else
                    if( (all(~ismember([CurPoint-1 CurPoint+2],BadTimePoints)) ) && ( CurPoint+2 <= size(CTC2DBig,2) ))
                        CTC2DBig(:,CurPoint)=mean(CTC2DBig(:,[CurPoint-1 CurPoint+2]),2);
                    else% Non of the if's occured
%                         error([num2str(CurPoint) ' -E- None of the above conditions are met!']);
                        AddToLog(WorkingP,'c_3badpnts',['Bad time points compensation: None of the conditions are met for ' num2str(BadTimePoints)],[],3);
                    end
                end
            end
        end
    end
end
%% Pad the CTC data with the T1 maps which are additional to main
SampleTsNoBefore=[(0:(NumVols-1))*TimeBetweenDCEVolsFinal];
% Make sure we actually have more than 1 T1 base map
if (Num_Of_T1_Maps>1)
    SampleTsNoBefore=[(0:(NumVols-1))*TimeBetweenDCEVolsFinal Additional_T1_Maps_Time_Diff_Sets(Additional_T1_Maps_Time_Diff_Sets>0)'*86400];
    %     kCoeffs=NaN(1,Num_Of_T1_Maps);
    FinalT1s=NaN(size(CAdditional_T1_Maps));
    %     for t=1:Num_Of_T1_Maps
    %         TmpT1=squeeze(CAdditional_T1_Maps(:,:,:,t));
    %         RefVolX=RefVol & TmpT1>100 & Baseline>100;
    %         S=sort(TmpT1(RefVolX));
    %         RefCurVal=S(max(1,floor(numel(S)*WhichPercent)));
    %         % RefCurVal=median(CT1(RefVolX));
    %         kCoeffs(t)=RefTrgVal/RefCurVal;
    % %         FinalT1s(:,:,:,t)=CAdditional_T1_Maps(:,:,:,t)*kCoeffs(t);
    % %         FinalT1s(:,:,:,t)=(Additional_T1_Maps(:,:,:,t).*PseudoB1.^2)*kCoeffs(t);
    %     end
    %     AddToLog(WorkingP,'c_20Norms',['kCoeffs ' num2str(kCoeffs)]);
    
    NotBaseIdxs=setdiff(1:Num_Of_T1_Maps,WhichSetForBase);
    RNotBaseIdxs=NotBaseIdxs*NaN;
    RNotBaseIdxs(NotBaseIdxs)=1:(Num_Of_T1_Maps-1);
    
    load(RelaxFNx,'TRsBySet');
    
    for qq=[WhichSetForBase SingleAngleIdxs]
        dir_to_look_name = dir([WorkingP 'Relaxometry' filesep 'Set_Num_' num2str(qq) filesep 'Core*.nii']);
        TmpFN=[WorkingP 'Relaxometry' filesep 'Set_Num_' num2str(qq) filesep dir_to_look_name(1).name];
        FullVol=loadniidata(TmpFN);
        FullVol4D(:,:,:,qq)=FullVol;
%         Tmp=DCEVolToT1(FA3D,DCEM0,FullVol);
        Tmp=DCEVolToT1TR(BaselineFA3Dx,DCEM0,FullVol,TRsBySet(qq));
        
        FinalT1s(:,:,:,qq)=Tmp;
        Additional_T1_Maps(:,:,:,qq)=Tmp;
        CAdditional_T1_Maps(:,:,:,qq)=Tmp;
    end
    
    DESPOT1Idxs=setdiff(1:Num_Of_T1_Maps,SingleAngleIdxs);
    for qq=DESPOT1Idxs
        dir_to_look_name = dir([WorkingP 'Relaxometry' filesep 'Set_Num_' num2str(qq) filesep 'Core*.nii']);
        Tmp=cell2mat({dir_to_look_name.name}');
        TmpFAs=str2num(Tmp(:,16:17));
        [AA,Chosen]=min(abs(TmpFAs-DCEFA));
        TmpFN=[WorkingP 'Relaxometry' filesep 'Set_Num_' num2str(qq) filesep dir_to_look_name(Chosen).name];
        FullVol=loadniidata(TmpFN);
        FullVol4D(:,:,:,qq)=FullVol;
%         Tmp=DCEVolToT1(FA3D,DCEM0,FullVol);
        Tmp=DCEVolToT1(BaselineFA3Dx,DCEM0,FullVol);
        
        FinalT1s(:,:,:,qq)=Tmp;
        Additional_T1_Maps(:,:,:,qq)=Tmp;
        CAdditional_T1_Maps(:,:,:,qq)=Tmp;
    end
    
    % Clean by PseudoB1
    %     FinalT1s=Additional_T1_Maps.*(repmat(PseudoB1.^2,[1 1 1 Num_Of_T1_Maps])).*kCoeff;
    %     FinalT1s=CAdditional_T1_Maps.*kCoeff;
    Raw2Nii(Additional_T1_Maps,[WorkingP 'UncleanedT1Maps.nii'],'float32', MeanFN);
    Raw2Nii(CAdditional_T1_Maps,[WorkingP 'CleanedT1Maps.nii'],'float32', MeanFN);
    Raw2Nii(FinalT1s,[WorkingP 'FinalT1Maps.nii'],'float32', MeanFN);
    %     FinalT1s=CAdditional_T1_Maps.*kCoeff;
    
    %     figure;montage(mritransform(repmat(FinalT1s(GoodRows,GoodCols,MidSli,:),[1 1 3 1])/3000));
    %     DiffT1=FinalT1s-repmat(FinalT1,[1 1 1 Num_Of_T1_Maps]);
    %     figure;montage(mritransform(repmat(500+DiffT1(GoodRows,GoodCols,MidSli,:),[1 1 3 1])/1000))
    
    %     Additional_T1_Inverse = 1 ./ (Additional_T1_Maps(:,:,:,NotBaseIdxs)*kCoeff);
    Additional_T1_Inverse = 1 ./ (FinalT1s(:,:,:,NotBaseIdxs));
    
    % Get the 2D additional maps
    AdditionalT1_Inverse_2D = Reshape4d22d(Additional_T1_Inverse,DBrainMask);
    
    % Calculate the base rate
    R10=1./FinalT1;
    % R0 = 1 / T0
    % Mask the unwanted voxels
    R10Clmn=R10(DBrainMask);
    
    % C(t) = 1/T1(t) - 1/(T1base)
    % Get the 2D additional maps
    Additional_CTC_2D=repPlus(AdditionalT1_Inverse_2D,-R10Clmn);
    %     for qq=SingleAngleIdxs
    %         dir_to_look_name = dir([WorkingP 'Relaxometry' filesep 'Set_Num_' num2str(qq) filesep 'Core*.nii']);
    %         TmpFN=[WorkingP 'Relaxometry' filesep 'Set_Num_' num2str(qq) filesep dir_to_look_name(1).name];
    %         FullVol=loadniidata(TmpFN);
    %         Tmp=DCEVolToT1(FA3D,DCEM0,FullVol);
    %         Tmp2=1./Tmp-1./FinalT1;
    %         Additional_CTC_2D(:,RNotBaseIdxs(qq))=Tmp2(DBrainMask);
    %
    % %         CAdditional_T1_Maps(:,RNotBaseIdxs(qq))=Tmp2(DBrainMask);
    % %         OAdditional_T1_Maps(:,RNotBaseIdxs(qq))=Tmp2(DBrainMask);
    %     end
    Additional_CTC_2D(imag(Additional_CTC_2D(:))~=0)=0;
    CTC2DBig(imag(CTC2DBig(:))~=0)=0;
    
    FigShowingCTCUsingDCEMainAndAdditional=false;
    if(FigShowingCTCUsingDCEMainAndAdditional)
        figure(9292);clf;
        Tmp=DBrainMask*0;
        Mn=-0.0005;
        Mx=0.0005;
        for qq=1:size(Additional_CTC_2D,2)
            Tmp(DBrainMask)=Additional_CTC_2D(:,qq);
            gsubplot(size(Additional_CTC_2D,2)+1,qq);
            imagesc(mritransform(Tmp(:,:,MidSli)),[Mn Mx])
            if(ismember(NotBaseIdxs(qq),SingleAngleIdxs))
                title('Single');
            end
        end
        gsubplot(size(Additional_CTC_2D,2)+1,size(Additional_CTC_2D,2)+1);
%         Tmp(DBrainMask)=CTC2DBig(:,end);
        Tmp(ImagB3D)=CTC2D(:,NumVols);
        imagesc(mritransform(Tmp(:,:,MidSli)),[Mn Mx])
        title('Last of DCEMain');
    end
    
    % Get the 2D additional maps
    
    %     num_before             = T1_Maps_Time_Vector(T1_Maps_Time_Vector < Main_Vector_Time(1));
    Additional_before_main = Additional_CTC_2D(:,intersect(1:(num_before-1),GoodAnglesF));
    Additional_after_main  = Additional_CTC_2D(:,intersect(num_before :end,GoodAnglesF));
    
%     CTC2DBig = [Additional_before_main CTC2DBig Additional_after_main];
    CTC2DBig = [CTC2DBig Additional_after_main];
end

T1All4D=Reshape2DCto4D(mat2cell(1./(CTC2DBig+repmat(1./FinalT1(DBrainMask>0),[1 size(CTC2DBig,2)])),size(CTC2DBig,1),ones(1,size(CTC2DBig,2))),DBrainMask>0);
Raw2Nii(T1All4D,[WorkingP 'T1All4D.nii'],'float32', MeanFN);
clear T1All4D
% % Reshape the 4d image back to 2d
% CTC2DBig = Reshape4d22d(CTC4DBig,DBrainMask);


%%

% ASK GILAD - what is this mask exactly and why did he use "imag"?
% ANSWER       - Matlab gives an imaginary number for certain calculations (such as log(-5)).This is non-physiological.
%                He takes a mask not including NANs and Imaginary numbers during the entire time of the test ( any("",2) ).
MskCTCGood=~any(imag(CTC2DBig)~=0,2) | any(isnan(CTC2DBig),2);
CTC2DBigGood=CTC2DBig(MskCTCGood(:),1:UnderSampling:end);
% save([WorkingP 'CTCMsk' USStr '.mat'],'CTC2DBigGood','MskCTCGood');

F1=find(DBrainMask);
F2=find(Msk);
B3=ismember(F1,F2);
CTC2DA=CTC2DBig(B3(:),1:UnderSampling:end);

clear CTC2DBig

% ASK GILAD - Ask him to explain all the masks till the AIF_Parker section.
ImagB=any(imag(CTC2DA)~=0,2) | any(isnan(CTC2DA),2);
CTC2DB=CTC2DA(~ImagB,:);
% MinusB=any(CTC2DB<-0.0001,2);
% CTC2D=CTC2DB(~MinusB,:);
CTC2D=CTC2DB;
if(isempty(CTC2D))
    error('Empty CTC2D');
end

Msk2=Msk;
Msk2(Msk2)=~ImagB;
% Msk2(Msk2)=~MinusB;
% Msk2(Msk2)=Baseline(Msk2)>10;
% DCE2D2=Reshape4d22d(DCE4D,Msk2);

ImagB3D=Msk;
ImagB3D(Msk)=~ImagB;

% MinusB3D=ImagB3D;
% MinusB3D(MinusB3D)=~MinusB;

% DCE2D=Reshape4d22d(DCE4D,Msk2);
% Baseline2D=Baseline(Msk2);
% RDCE2D=repMulti(DCE2D,1./Baseline2D);
% Msk2F=find(Msk2(Msk));size(Msk2F,1)

% doesn't use FA3D
% T12D=CT1(Msk2);
% R102D=R10(Msk2);
% Tt2D=TtFunc(repmat(T12D,1,nVols), RDCE2D, DCETR,DCEFA);
% CTC2D=1./Tt2D-repmat(R102D,1,nVols);

% ASK GILAD - What are the following values and how did he determine it? (they are a bit different than what I see in the article)
% ASNWER - Parker, Geoff JM, et al. "Experimentally?derived functional form for a population?averaged high?temporal?resolution arterial input function for dynamic contrast?enhanced MRI." Magnetic Resonance in Medicine 56.5 (2006): 993-1000.
A1=0.809;A2=0.330;T1=0.17046;T2=0.364;
sig1=0.055;sig2=0.134;alphaa=1.064;beta=0.166;
s=37.772;tau=0.482;
% Before there were really slighltly different numbers, from other place?
% A~mM.min, T~min; sig~min, alpha~mM, beta,s~1/min, tau~min
% Time stamp (in minutes) for every temporal point
ts=((1:nVols)-3)./Min2Timestep;

% Population average c(t) according to Parker's article.
C=AIF_Parker(ts,A1,sig1,T1,A2,sig2,T2,alphaa,beta,s,tau)/1000;
C(1:4)=0;
% ASK GILAD - why isnt he interpolating eventually?
% ANSWER - He is interpolating it, but later on (AIF finding + PK extraction).
TimeBetweenDCEVolsFinal=TimeBetweenDCEVolsFinal*UnderSampling;
% Time between DCE vols in seconds is used for determining interpolation
InterpolationFactor=ceil(TimeBetweenDCEVolsFinal);
% tic;
% SICTC2D=SmoothAndInterpolate(CTC2D,InterpolationFactor,BolusStart);
% t=toc;
% disp(['Smoothing took ' num2str(t) 's']);
%%
CTC4D=Reshape2DCto4D(mat2cell(CTC2D,size(CTC2D,1),ones(1,size(CTC2D,2))),Msk2);
Raw2Nii(CTC4D,[WorkingP 'CTC4D.nii'],'float32', MeanFN);
%%
clear FullVol4D T12D IF T1FromBaseline Tmp
NumVols=size(CTC2D,2);
DCE4DX=max(0,CTC4D);
[MaxVal, PeakTime]=max(DCE4DX,[],4);
PeakTime(~FBrainMask)=NaN;
GoodPeak=PeakTime>(BolusStart-5) & PeakTime<(BolusStart+15);
[I J K]=ind2sub(size(GoodPeak),find(GoodPeak));
NearVals=[DCE4DX(sub2ind(size(DCE4D),I,J,K,max(1,PeakTime(GoodPeak)-1))) DCE4DX(sub2ind(size(DCE4D),I,J,K,min(NumVols,PeakTime(GoodPeak)+1)))];
MaxVal1D=MaxVal(GoodPeak);
PeakTime1D=PeakTime(GoodPeak);
[MaxNearVal, WhichDirection]=max(NearVals,[],2);
RelativeDist=MaxNearVal./(MaxVal1D+MaxNearVal);
ApproximatePeakTime1D=PeakTime(GoodPeak)+((WhichDirection-1)*2-1).*RelativeDist;
ApproximatePeakTime=PeakTime;
ApproximatePeakTime(GoodPeak)=ApproximatePeakTime1D;
ApproximatePeakTime(~GoodPeak)=NaN;
% MoreExactPeakTime(~GoodPeak)=NaN;
clear DCE4DX CTC4D T12D IF T1FromBaseline Tmp
clear BaselineB1Fld Baseline BaselineFA3D CurT1 B1Nok FA3DNoB1Nok CurVol
clear Additional_T1_Inverse Additional_T1_Maps CAdditional_T1_Maps
%% Clearing vars
clear DCE4D pnk mixture Tt2D DCE2D RDCE2D CTC2DBrainMsk
clear FA4D AD3DAM ARelax2Step ARelaxFAsearch ARelaxWndFA E1F TmpVar CTC2DB CTC2DA
%% Saving .mat file for that stage
save(CTCFN);
disp('DCET1_CTC finished');