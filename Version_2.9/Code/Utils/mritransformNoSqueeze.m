function Out=mritransformNoSqueeze(I)
I=permute(I,[2 1 3:numel(size(I))]);
Out=(flipdim(I,1));
