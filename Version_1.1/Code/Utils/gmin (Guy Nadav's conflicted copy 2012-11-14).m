function [Out Idxs]=gmin(in,dims,opts,Msk)
if(nargin>3)
    in(~Msk)=Inf;
end
if(~exist('dims','var'))
    [Out a]=min(in(:));
    if(nargout>1)
        [Idxs{1:ndims(in)}]=ind2sub(size(in),a);
        Idxs=[Idxs{:}];
    end
    return;
end
Out=min(in,[],dims(1));
for i=2:length(dims)
    Out=min(Out,[],dims(i));
end
if(~exist('opts','var'))
    return;
end
if(ismember(1,opts))
    Out=squeeze(Out);
end