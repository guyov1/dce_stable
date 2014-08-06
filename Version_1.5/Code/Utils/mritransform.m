function Out=mritransform(I)
I=permute(I,[2 1 3:numel(size(I))]);
Out=squeeze(flipdim(I,1));
