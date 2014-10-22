% function Out=SPGRfM(T1s,PDs,FAs,TR)
function Out=SPGRfM(T1s,PDs,FAs,TR)
nFAs=numel(FAs);
n=numel(T1s);
TRM=repmat(TR,[n 1]);
T1M=repmat(T1s',[1 nFAs]);
PDM=repmat(PDs',[1 nFAs]);
if (numel(FAs) == 1)
    
    if ( size(TRM,1) == 1 )
        TRM = transpose(TRM);
    end
    if ( size(T1M,1) == 1 )
        T1M = transpose(T1M);
    end
    if ( size(PDM,1) == 1 )
        PDM = transpose(PDM);
    end
    
    E1=exp(-TRM./T1M);
    
    SFAM=repmat(sind(FAs),[1 n])';
    CFAM=repmat(cosd(FAs),[1 n])';
    
    
else
    E1=exp(-TRM./T1M);
    SFAM=repmat(sind(FAs),[n 1]);
    CFAM=repmat(cosd(FAs),[n 1]);

end


Out=PDM.*(1-E1).*SFAM./(1-E1.*CFAM);