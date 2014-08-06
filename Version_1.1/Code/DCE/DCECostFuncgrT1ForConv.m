% function Out=DCECostFunc(AIF,Ktrans,Ve,Vp,TimePoints,ConvIdxMTriB,TR,CosFA,R10,r1,Ws,RData)
function Out=DCECostFuncgrT1ForConv(AIF,Kep,TimePoints,ConvIdxMTriB,TriB)

% Number of possible Keps 
n=numel(Kep);

% Number of time points
nTimePoints=numel(TimePoints);

% CosFA=cosd(FA);

% Interval between 2 time points
dT=TimePoints(2)-TimePoints(1);

% Avoid the same Kep
[U, QQQ, IB]=unique(Kep);

% Number of unique Keps
nUKep=numel(U);

% Create a zero matrix in the size of the convolution matrix
EKepM=zeros(nTimePoints);

% TriB=ConvIdxM>0;

% Zero matrix for the convolution result 
CAIF=zeros(nTimePoints,nUKep);

% For each possible kep
for i=1:nUKep
    
    %     - Kep*t
    %   e
    % 
    % For all possible t's -> [ exp(-kt1)  exp(-kt2)  exp(-kt3)  ....exp(-ktN)  ]
    EKepV=exp( -U(i) * TimePoints);
    
    %     EKepM(TriB)=EKepV(ConvIdxM(TriB));
    
    % Fill the convolution matrix with the exp(-kep*t) values
    %
    %   [ exp(-kt1)      0                         0                       ....    0                   ]
    %   [ exp(-kt2)      exp(-kt1)        0                       ....    0                   ]
     %  [                               ............                                                               ]  
    %   [ exp(-ktN-1)  exp(-ktN-2)    exp(-ktN-3)  ....    0                   ]
    %   [ exp(-ktN)     exp(-ktN-1)    exp(-ktN-2)  ....    exp(-kt1)  ]
    %
    EKepM(TriB)=EKepV(ConvIdxMTriB);
    
    % Holds the result of the convolution for the specific Kep
    %        - Kep*t
    %   e                       CONV   AIF(t)
    %
    CAIF(:,i)=EKepM*AIF;
    
end

% Create an index matrix for the convoultion result matrix
%
% [ 1    2         ... N ]
% [N    N+1   ... 2N]
% [2N 2N+1  ... 3N]
% [          ...            ]
%
IdxMat=repmat(IB-1,1,nTimePoints) * nTimePoints+repmat(1:nTimePoints,n,1);

% get C
% ASK GILAD - What is the meaning of the following condition? Seems like
%                           both are the same...
% ASNWER - The first one has transpose on it. Used when we work on only one voxel.
% Multiply by dT to get the real time stamps (and not indices 1, 2, 3 ...)
if(numel(IB)==1)
    Out=CAIF(IdxMat)'*dT;
else
    Out=CAIF(IdxMat)*dT;
end