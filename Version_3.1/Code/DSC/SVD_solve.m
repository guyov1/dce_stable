function Rt = SVD_solve_voxel(AIF,Ct,deltaT,method)

% This function solve (Ab=c) using SVD for b, (as shown in the paper of
% Ostergard 1996), for a single voxel. 
%   A is matrix built from AIF, 
%   Ct=c=concentration curve at a voxel, method='constant' or 'linear' (assumption on the
%       difference between 2 time points),
%   b=Rt (the desired solution).
% 

AIF=AIF(:);
AIFmat=AIF2mat(AIF,deltaT,method);

if length(AIF)~=length(Ct)
    error('AIF and concentration curve are not the same length');
end

[U,W,V] = svd(AIFmat); % A=U*W*V'
W_inv=inv(W);
W_inv_after_proc=process_W_SVD(W_inv);  % here we can "regularize" the problem by eliminating elements from diagonal

Rt= v*W_inv_after_proc*U'*Ct(:); %  b= VWU'c
