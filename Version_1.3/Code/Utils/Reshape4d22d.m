function Out2d=Reshape4d22d(M4d,B3d)

% Time dimension
N=size(M4d,4);

% Get row, column and slice index of all relevant indices determined by the mask
[I J K]=ind2sub(size(B3d),find(B3d));

% Arrange the indices in the following way:
%
%                         x1 y1 z1     1
%                         x1 y1 z1     2
%                         .......................
%                         x1 y1 z1     50 (time frames)
%                         x2 y2 z2     1
%                         x2 y2 z2     2
%                         .......................
%                         x2 y2 z2     50 (time frames)
%                        ........................
%                        ........................
%                         xN yN zN     1
%                         xN yN zN     2
%                         .......................
%                         xN yN zN     50 (time frames)

Idxs=[kron([I J K],ones(N,1)) repmat((1:N)',numel(I),1)];

% Get the linear index of all voxels during time (from 4D)
NIdxs=sub2ind(size(M4d),Idxs(:,1),Idxs(:,2),Idxs(:,3),Idxs(:,4));

% Reshape the voxels in the 4D map (according to the needed indices) to 2D
Out2d=reshape(M4d(NIdxs)',N,numel(I))';