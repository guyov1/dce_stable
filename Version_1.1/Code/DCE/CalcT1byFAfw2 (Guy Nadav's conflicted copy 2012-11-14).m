function [T1s PDs RMS]=CalcT1byFAfw2(Data,FAs,TRs)
% Data - nFAs x nSamples
% FAs - Row vector
T1Base=1200;
T1Max=40000;
TR=TRs(1);
SinFAs=sind(FAs)';
CosFAs=cosd(FAs)';
n=size(Data,2);

C=Data.*((repmat(CosFAs,1,n)));
X=zeros(2,n);
CurM=[C(:,1) SinFAs];
% RMS=zeros(n,1);
W=1./(1-exp(-TR/T1Base)*CosFAs);

Wc=repmat(W,1,2);

for i=1:n
    CurM(:,1)=C(:,i);
    XTmp=(CurM.*Wc)\(Data(:,i).*W);
%     TTmp=-TR./log(XTmp(1));
%     W2=1./(1-exp(-TR/TTmp)*CosFAs);
    W2=1./(1-XTmp(1)*CosFAs);
    Wc2=repmat(W2,1,2);
    X(:,i)=(CurM.*Wc2)\(Data(:,i).*W2);

%     RMS(i)=norm( Data(:,i)-CurM*X(:,i) );
end
X=X';
T1s=-TR./log(X(:,1));
PDs=X(:,2)./(1-X(:,1));
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
Diffs=SPGRf(T1s,PDs,FAs,TR)'-Data;
RMS=sqrt(mean(Diffs.^2,1));