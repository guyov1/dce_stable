function Out2d=Reshape4d22d(M4d,B3d)
N=size(M4d,4);
[I J K]=ind2sub(size(B3d),find(B3d));
Idxs=[kron([I J K],ones(N,1)) repmat((1:N)',numel(I),1)];
NIdxs=sub2ind(size(M4d),Idxs(:,1),Idxs(:,2),Idxs(:,3),Idxs(:,4));
Out2d=reshape(M4d(NIdxs)',N,numel(I))';