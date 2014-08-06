% function Out=SPGRf(T1s,PDs,FAs,TR)
% Calculate the signal out of the estimated parameters
function Out=SPGRf(T1s,PDs,FAs,TR)



nFAs=numel(FAs);
n=numel(T1s);


E1=repmat(exp(-TR./T1s),1,nFAs);

%               -TR
%                --
%                T1
%     Pd * (1 - e  ) * sin(a)
%     -----------------
%                   -TR
%                    --
%                    T1
%       1  - cos(a)*e

Out=repmat(PDs,1,nFAs).*(1-E1).*repmat(sind(FAs),n,1)./(1-E1.*repmat(cosd(FAs),n,1));