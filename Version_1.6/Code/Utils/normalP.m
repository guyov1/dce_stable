function [out minval maxval]=normalP(in,pl, ph)
if(~exist('pl','var'))
    pl=[];
end
if(~exist('ph','var'))
    ph=[];
end
if(isempty(pl))
    pl=0.1;
end
if(isempty(ph))
    ph=0.9;
end
S=sort(in(isfinite(in) & in~=0));
N=numel(S);
if(N<=1)
    minval=0;
    maxval=1;
else
    minval=S(max(1,floor(N*pl)));
    maxval=S(ceil(N*ph));
end
range=maxval-minval;
out=(in-minval)./range;