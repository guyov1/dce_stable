function [T1s PDs RMS]=CalcT1byFAfw2_Single_Angle(Data,FAs,TRs,M0_Map)
% Data - nFAs x nSamples
% FAs  - Row vector
% TRs  - Row vector

% ASK GILAD - why did he pick those values?
T1Base = 1200;
T1Max  = 40000;

% Get the first TR
TR=TRs(1);

% Get the sine and the cosine
SinFAs=sind(FAs)';
CosFAs=cosd(FAs)';

% Number of data points (voxels over time)
n=size(Data,2);

% Data multiplied by cos(flip angle)
S_Cos  = Data .* repmat(CosFAs,1,n);
% M0 multiplied by sin(flip angle)
M0_Sin = M0_Map' .* repmat(SinFAs,1,n);

%                   -TR
%           ---------------------
% T1 = 
%           (   S - M0*sin(a)     )
%           (  -----------------  )
%        ln (                     )
%           (  S*cos(a)-M0*sin(a) )

X   =  (Data - M0_Sin) ./ (S_Cos - M0_Sin);
T1s = -TR ./ log(X);

PDs = M0_Map';

F=find(imag(T1s)~=0 | imag(PDs)~=0);
T1s=max(100,min(T1s,T1Max));
T1s(F)=0;
PDs(F)=0;


% Calculate the signal out of the estimated parameters and measure the error
Diffs=SPGRf(T1s,PDs,FAs,TR)-Data;
RMS=sqrt(mean(Diffs.^2,1));


