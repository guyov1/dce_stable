function Out=Reshape2DCto4D(in,Msk)

% If it is not a cell, make it a cell
if(~iscell(in))
    in={in};
end

Out=zeros([size(Msk) numel(in)]);
Tmp=Msk*0;

for i=1:numel(in)
    Tmp(Msk)=in{i};
    if(ndims(Msk)==3)
        Out(:,:,:,i)=Tmp;
    else
        Out(:,:,i)=Tmp;
    end
end
Out=squeeze(Out);