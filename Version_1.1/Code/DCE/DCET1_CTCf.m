% The SPGR equation:
% S(t)=M0*sin(FA)*(1-E1(t))/(1-E1(t)cos(FA))
% where E1=exp(-TR/T1(t))
% So, given M0 from relaxometry (and TR, FA, S), we can get T1(t)
% The contrast agent concentration equation
% 1/T1(t)=1/(T1base)+C(t)
% So given T1base, we can get from T1(t) the Concentration Time Curve C(t)
function DCET1_CTCf(DCEInfo,WorkingP,DoN3,DoGlobal,DoDeviations,CalcForce,Options)
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
disp('DCE_CTC Starting');

RWorkingP=[WorkingP 'Relaxometry/'];

FMaskFN=[RWorkingP 'FBrainMsk.nii'];

T1Res{1,1,1}=[RWorkingP 'T13DOFA.nii'];
T1Res{1,2,1}=[RWorkingP 'T13DNFA.nii'];

% ASK GILAD - Didn't he run the b1 cleaning in DCET1_RelaxForSubjectf?
if(DoN3)
    T1Res{1,1,2}=Iterative_n3_b1_clean_mask_FN(3,T1Res{1,1,1},FMaskFN,false,[]);
    if(DoDeviations)
        T1Res{1,2,2}=Iterative_n3_b1_clean_mask_FN(3,T1Res{1,2,1},FMaskFN,false,[]);
    end
end
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

%% Setting threshold values
% ASK GILAD - how did he determine those values?
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

%% M0
load(PrepareFN,'nVols','Min2Timestep','BrainMask');

% Get FA and TR
DCETR=DCEInfo.RepetitionTime;
DCEFA=DCEInfo.FlipAngle;

% Data saved at the end of DCET1_RelaxForSubjectf.m
MatFN=[RWorkingP 'NFARes.mat'];
load(MatFN);

% ASK GILAD - I don't understand the B1 calculations
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
% ASK GILAD - how does he calculate FA and what is the purpose of the term
%             CT1*0?
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
DCERelaxP=[WorkingP 'Relaxometry/'];
FMaskFN=[DCERelaxP 'FBrainMsk.nii'];
FBrainMask=loadniidata(FMaskFN)>0;
% ASK GILAD - what is the purpose of the following 2 lines?
se=strel('disk',30,8);
DBrainMask=imdilate(FBrainMask,se);

% Compute two dimenstional C(t)
% Note - The following is a heavy computation and can freeze the computer for a minute
CTC2DBig=DCEtoCTCf(DBrainMask,DCE4D(:,:,:,1:(end-Options.nVolsToRemoveFromEnd)),CT1,FA3D,DCETR, Baseline);

% ASK GILAD - what is this mask exactly and why did he use "imag"?
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

%% ASK GILAD - what are the following values and how did he determine it? (they are a bit different than what I see in the article)
A1=0.809;A2=0.330;T1=0.17046;T2=0.365;
sig1=0.0563;sig2=0.132;alpha=1.050;beta=0.1685;
s=38.078;tau=0.483;
% A~mM.min, T~min; sig~min, alpha~mM, beta,s~1/min, tau~min
% Time stamp (in minutes) for every temporal point
ts=((1:nVols)-3)./Min2Timestep;

% Population average c(t) according to Parker's article.
C=AIF_Parker(ts,A1,sig1,T1,A2,sig2,T2,alpha,beta,s,tau)/1000;
C(1:4)=0;
%% ASK GILAD - why isnt he interpolating eventually?
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