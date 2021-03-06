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

% Starting CTC calculation message
disp('DCE_CTC Starting');

% Get the T1s data from the first set calculated
RWorkingP=[WorkingP 'Relaxometry' filesep 'Set_Num_1' filesep];
FMaskFN=[RWorkingP 'FBrainMsk.nii'];

% Get the T1 base map (the second one is FA corrected - if exists)
T1Res{1,1,1}=[RWorkingP 'T13DOFA.nii'];
T1Res{1,2,1}=[RWorkingP 'T13DNFA.nii'];

% ASK GILAD - Didn't he run the b1 cleaning in DCET1_RelaxForSubjectf?
% ASNWER  - He used it just to get the file name (T1Res{1,1,2}).
%                      If it was already calculated, it won't be calculated again(using "false" option).
if(DoN3)
    T1Res{1,1,2}=Iterative_n3_b1_clean_mask_FN(3,T1Res{1,1,1},FMaskFN,false,[]);
    if(DoDeviations)
        T1Res{1,2,2}=Iterative_n3_b1_clean_mask_FN(3,T1Res{1,2,1},FMaskFN,false,[]);
    end
end

% Get the N3 corrected T1 maps
T1Res{1,1,3}=[RWorkingP 'T13DOFA_N3k.nii'];
T1Res{1,2,3}=[RWorkingP 'T13DNFA_N3k.nii'];

RelaxFN=T1Res{1,DoDeviations+1,DoN3+DoGlobal+1};

% Load .nifti data of the T1 map
CT1=loadniidata(RelaxFN);

% Get the mean DCE map
MeanFN=[WorkingP 'DCEMean.nii'];

% ASK GILAD - what is the purpose of the following?
if(DoN3)
    PseudoB1=sqrt(double(loadniidata(T1Res{1,DoDeviations+1,2}))./double(loadniidata(T1Res{1,DoDeviations+1,1})));
    B1FN=[RWorkingP 'N3B1.nii'];
    Raw2Nii(PseudoB1,[RWorkingP 'N3B1.nii'],'float32', MeanFN);
end

%% Load all previous stages
PrepareFN=[WorkingP 'AfterPrepare4D.mat'];
load(PrepareFN,'Baseline','Msk','DCE4D','TimeBetweenDCEVols','BolusStart');

%% Update 4D map according to additional DESPOT1 data we have

% DCERelaxP=[WorkingP 'Relaxometry' filesep 'Set_Num_1' filesep];

% Get all T1 maps we have
All_T1_Maps = dir( [WorkingP 'Relaxometry' filesep 'Set*']);
% Get number of T1 maps
Num_Of_T1_Maps = numel(All_T1_Maps);
% If we have more than 1, use it and create time vector
if (Num_Of_T1_Maps>1)

    % Initiate time vector
    T1_Maps_Time_Vector = zeros(1,Num_Of_T1_Maps-1);
    % Initiate date vector
    T1_Maps_Time_Date      = zeros(1,Num_Of_T1_Maps-1);
    % Initiate T1 additional maps
    Additional_T1_Maps = zeros([size(CT1) Num_Of_T1_Maps-1]);

    % Ignore the first map which is used for T1-base calculation
    for map_idx = 2:Num_Of_T1_Maps
        dir_to_look = All_T1_Maps(map_idx);
        dir_to_look_name = [WorkingP 'Relaxometry' filesep dir_to_look.name filesep 'Middle_Time_Stamp.mat']; %#ok<AGROW>
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
        T1_Maps_Time_Vector(map_idx-1) = time_seconds;
        T1_Maps_Time_Date(map_idx-1)      = str2double(Mid_Time_Stamp_Date);

        % Update T1 nifti maps
        T1_Path = [WorkingP 'Relaxometry' filesep dir_to_look.name filesep];
        if (DoN3+DoGlobal+1 == 3)
            if (DoDeviations+1 == 2)

                T1_Map_Path=[T1_Path 'T13DNFA_N3k.nii'];
            else
                T1_Map_Path=[T1_Path 'T13DOFA_N3k.nii'];
            end
        else
            if (DoDeviations+1 == 2)
                T1_Map_Path=[T1_Path 'T13DNFA.nii'];

            else
                T1_Map_Path=[T1_Path 'T13DOFA.nii'];
            end

        end

        Additional_T1_Maps(:,:,:,map_idx-1) = loadniidata(T1_Map_Path);

    end

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
    if (length(unique([T1_Maps_Time_Date Main_First_Time_Stamp_Date])) ~= 1)
        error(' Error! - The dates of the T1 maps are different! Probably did it around midnight!');
    end

    % Create the time vector for the main acquisition
    % I assumes the acquisitions are uniformly distributed from the first
    % main image time to the first image after the main time
    all_bigger_than_main_time = T1_Maps_Time_Vector(T1_Maps_Time_Vector>time_seconds_main);
    first_after_main_time = all_bigger_than_main_time(1);
    Main_Vector_Time = round(linspace(time_seconds_main,first_after_main_time,num_volumes+1));
    Main_Vector_Time = Main_Vector_Time(1:end-1);
end


%DCE4D

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

    % If the number of pixels in slice, exceeding T1Thresh is bigger than NaTThresh, mark the slice as NaN.
    if(nAT(i)>NaTThresh)
        CT1(:,:,i)=NaN;
    end

end

%% Calculate M0

load(PrepareFN,'nVols','Min2Timestep','BrainMask');

% Get FA and TR
DCETR=DCEInfo.RepetitionTime;
DCEFA=DCEInfo.FlipAngle;

% Data saved at the end of DCET1_RelaxForSubjectf.m
MatFN=[RWorkingP 'NFARes.mat'];
load(MatFN);

% ASK GILAD - I don't understand the B1 calculations
% ASNWER       - SHOULD READ THE ARTICLE AND CHECK THAT.
if(DoN3)
    B1FN=[RWorkingP 'N3B1.nii'];
    B1=loadniidata(B1FN);
else
    B1=CT1*0+1;
end

if(DoGlobal)
    B1=B1*sqrt(kCoeff(1+DoDeviations));
end

% Basic rate  - 1/T1base
R10=1./CT1;
% E1=exp(-TR/T1(t))
E1=exp(-DCETR.*R10);

% Flip Angle
% ASK GILAD - How does he calculate FA and what is the purpose of the term CT1*0?
% ANSWER       - He used the CT1*0 just to use the matrix dimensions of CT1.
%                          He multiplied by B1 to take into consideration the fact the angle is not uniform (because of B1).
FA3D=CT1*0+DCEFA;
FA3D=FA3D.*B1;

% Cos and Sin of Flip Angle
CosFA=cosd(FA3D);
SinFA=sind(FA3D);

% M0 =  ( S(t0) * (1-E1(t)cos(FA)) ) / sin(FA)*(1-E1(t))
% The base line is the mean of the first images until the bolus (DCET1_Prepare4Df.m)
DCEM0=(Baseline.*(1-E1.*CosFA))./((1-E1).*SinFA);
DCEM0FN=[WorkingP 'DCE_M0.nii'];
Raw2Nii(DCEM0,DCEM0FN,'float32',MeanFN);

% M0Clmn=DCEM0(Msk);
% %%
% DCE2D=Reshape4d22d(DCE4D,Msk);
%
% TmpVar=repMulti(DCE2D,1./M0Clmn);
% E1F=(repPlus(TmpVar,-SinFA(Msk)))./(repPlus(repMulti(TmpVar,CosFA(Msk)),-SinFA(Msk)));
% R1F=-log(E1F)./DCETR;
% R10Clmn=R10(Msk);
% CTC2DA=repPlus(R1F,-R10Clmn);

% Get the T1s data from the first set calculated
DCERelaxP=[WorkingP 'Relaxometry' filesep 'Set_Num_1' filesep];
FMaskFN=[DCERelaxP 'FBrainMsk.nii'];
% Load the brain mask
FBrainMask=loadniidata(FMaskFN)>0;

% ASK GILAD - what is the purpose of the following 2 lines?
% ANSWER       - He wanted to take a little bigger mask than the calculated one.
se=strel('disk',30,8);
DBrainMask=imdilate(FBrainMask,se);

% Compute two dimenstional C(t)
% Note - The following is a heavy computation and can freeze the computer for a minute
CTC2DBig=DCEtoCTCf(DBrainMask,DCE4D(:,:,:,1:(end-Options.nVolsToRemoveFromEnd)),CT1,FA3D,DCETR, Baseline);


%% Pad the CTC data with the T1 maps which are additional to main

% Make sure we actually have more than 1 T1 base map
if (Num_Of_T1_Maps>1)
    
    Additional_T1_Inverse = 1 ./ Additional_T1_Maps;
    % Get the 2D additional maps
    AdditionalT1_Inverse_2D = Reshape4d22d(Additional_T1_Inverse,DBrainMask);
    
    % Calculate the base rate
    R10=1./CT1;
    % R0 = 1 / T0
    % Mask the unwanted voxels
    R10Clmn=R10(DBrainMask);
    
    % C(t) = 1/T1(t) - 1/(T1base)
    % Get the 2D additional maps
    Additional_CTC_2D=repPlus(AdditionalT1_Inverse_2D,-R10Clmn);

    % Get the 2D additional maps
    %Additional_CTC_2D = Reshape4d22d(Additional_CTC,DBrainMask);

    % T1_Maps_Time_Vector;
    % Main_Vector_Time;

    % Create the complete time vector (used in FindPKBATgAIFMurase.m)
    Complete_Time_Vector   = sort([T1_Maps_Time_Vector Main_Vector_Time]);
    
    num_before             = T1_Maps_Time_Vector(T1_Maps_Time_Vector < Main_Vector_Time(1));
    Additional_before_main = Additional_CTC_2D(:,1:num_before);
    Additional_after_main  = Additional_CTC_2D(:,num_before + 1:end);

    CTC2DBig = [Additional_before_main CTC2DBig Additional_after_main];

end 

% % Reshape the 4d image back to 2d
% CTC2DBig = Reshape4d22d(CTC4DBig,DBrainMask);


%%

% ASK GILAD - what is this mask exactly and why did he use "imag"?
% ANSWER       - Matlab gives an imaginary number for certain calculations (such as log(-5)).This is non-physiological.
%                           He takes a mask not including and NANs and Imaginary numbers during the entire time of the test ( any("",2) �).
MskCTCGood=~any(imag(CTC2DBig)~=0,2) | any(isnan(CTC2DBig),2);
CTC2DBigGood=CTC2DBig(MskCTCGood(:),1:UnderSampling:end);
save([WorkingP 'CTCMsk' USStr '.mat'],'CTC2DBigGood','MskCTCGood');

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
% ASNWER -        He does not remember. Maybe he took it from a similar article.
A1=0.809;A2=0.330;T1=0.17046;T2=0.365;
sig1=0.0563;sig2=0.132;alpha=1.050;beta=0.1685;
s=38.078;tau=0.483;
% A~mM.min, T~min; sig~min, alpha~mM, beta,s~1/min, tau~min
% Time stamp (in minutes) for every temporal point
ts=((1:nVols)-3)./Min2Timestep;

% Population average c(t) according to Parker's article.
C=AIF_Parker(ts,A1,sig1,T1,A2,sig2,T2,alpha,beta,s,tau)/1000;
C(1:4)=0;
% ASK GILAD - why isnt he interpolating eventually?
% ANSWER - He is interpolating it, but later on (AIF finding + PK extraction).
TimeBetweenDCEVols=TimeBetweenDCEVols*UnderSampling;
% Time between DCE vols in seconds is used for determining interpolation
InterpolationFactor=ceil(TimeBetweenDCEVols);
% tic;
% SICTC2D=SmoothAndInterpolate(CTC2D,InterpolationFactor,BolusStart);
% t=toc;
% disp(['Smoothing took ' num2str(t) 's']);


%% Clearing vars
clear DCE4D pnk mixture Tt2D DCE2D RDCE2D CTC2DBrainMsk
clear FA4D AD3DAM ARelax2Step ARelaxFAsearch ARelaxWndFA E1F TmpVar CTC2DB CTC2DA

%% Saving .mat file for that stage
save(CTCFN);
disp('DCET1_CTC finished');