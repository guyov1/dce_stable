% function Out=SPGRf(T1s,PDs,FAs,TR)
function Out=SPGRf(T1s,PDs,FAs,TR)
nFAs=numel(FAs);
n=numel(T1s);
E1=repmat(exp(-TR./T1s),1,nFAs);
Out=repmat(PDs,1,nFAs).*(1-E1).*repmat(sind(FAs),n,1)./(1-E1.*repmat(cosd(FAs),n,1));