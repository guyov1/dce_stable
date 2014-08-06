function [T1s PDs RMS]=CalcT1byFAfw2(Data,FAs,TRs)
% Data - nFAs x nSamples
% FAs - Row vector
% TRs - Row vector

% ASK GILAD - why did he pick those values?
T1Base=1200;
T1Max=40000;

% Get the first TR
TR=TRs(1);

% Get the sine and the cosine
SinFAs=sind(FAs)';
CosFAs=cosd(FAs)';

% Number of data points (voxels over time)
n=size(Data,2);

% Data multiplied by cos(flip angle)
C = Data .* repmat(CosFAs,1,n);
X = zeros(2,n);

% First column of data*cos(a) and sin(a)
CurM=[C(:,1) SinFAs];

% RMS=zeros(n,1);

%            1
%     -----------------
%                   -TR
%                    --
%                    T1
%       1  - cos(a)e

W = 1 ./ ( 1 - exp(-TR/T1Base) * CosFAs );

% Duplicate W to 2 columns
Wc=repmat(W,1,2);

% ASK GILAD - How does the following calculate T1??
% Go over each voxel
for i=1:n
    
    CurM(:,1)=C(:,i);
    
    XTmp = (CurM .* Wc) \ (Data(:,i) .* W);
    
%     TTmp=-TR./log(XTmp(1));
%     W2=1./(1-exp(-TR/TTmp)*CosFAs);

    W2  = 1 ./ (1-XTmp(1)*CosFAs);
    
    Wc2 = repmat(W2,1,2);
    
    X(:,i)= (CurM .* Wc2) \ (Data(:,i).*W2) ;

%     RMS(i)=norm( Data(:,i)-CurM*X(:,i) );
end


X=X';
T1s = -TR ./ log(X(:,1));
% ASK GILAD - What is PD?
PDs = X(:,2) ./ (1-X(:,1));
F=find(imag(T1s)~=0 | imag(PDs)~=0);
T1s=max(100,min(T1s,T1Max));
T1s(F)=0;
PDs(F)=0;


% XL=X*0;
% XL(:,1)=exp(-TR./T1s);
% XL(:,2)=PDs.*(1-X(:,1));
% XL=XL';
% RMSL=zeros(n,1);
% for i=1:n
%     CurM(:,1)=C(:,i);
%     RMSL(i)=norm( Data(:,i)-CurM*XL(:,i) );
% end

% Calculate the signal out of the estimated parameters and measure the error
Diffs=SPGRf(T1s,PDs,FAs,TR)'-Data;
RMS=sqrt(mean(Diffs.^2,1));