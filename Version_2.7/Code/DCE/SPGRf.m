% function Out=SPGRf(T1s,PDs,FAs,TR)
% Calculate the signal out of the estimated parameters
function Out=SPGRf(T1s,PDs,FAs,TR)



nFAs = numel(FAs);
nT1s = numel(T1s);


E1 = repmat(exp(-TR./T1s),1,nFAs);

%               -TR
%                --
%                T1
%     Pd * (1 - e  ) * sin(a)
%     -----------------
%                   -TR
%                    --
%                    T1
%       1  - cos(a)*e

if ( nFAs == 1 )
    Out = repmat(PDs,1,nFAs).*(1-E1).*repmat(sind(FAs),1,nT1s)./(1-E1.*repmat(cosd(FAs),1,nT1s));
else
    Out = repmat(PDs,1,nFAs).*(1-E1).*repmat(sind(FAs),nT1s,1)./(1-E1.*repmat(cosd(FAs),nT1s,1));
end
