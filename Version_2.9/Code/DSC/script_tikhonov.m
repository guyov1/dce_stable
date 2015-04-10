
% load('AIFmat.mat');
% load('Ct.mat');

% We use the Tikhonov regularization method for the problem:
% min( ||Ax-b||^2+lambda^2*||Lx||^2)
L1=zeros(size(AIFmat,1)-1,size(AIFmat,2));
zero_vec=zeros(size(Ct));

L1_main_diagonal=(-1)*ones(size(L1,2),1);
L1_second_diagonal=ones(size(L1,2)-1,1);

L1a=diag(L1_main_diagonal,0);
L1b=diag(L1_second_diagonal,1);

L1=L1a(1:end-1,:)+L1b(1:end-1,:);

% generalized SVD:
[U,V,X,C,S] = gsvd(AIFmat,L1);
mu=diag(S);
s=diag(C);
sigma=gsvd(AIFmat,L1);
% sigma=sqrt(diag(C'*C)./diag(S'*S));
sm=[reshape(sigma(1:end-1),length(sigma)-1,1) mu(:)];

% plot the L-curve:
% [reg_corner,rho,eta,reg_param]=l_curve(U,s,Ct,'Tikh',L1,V);
[reg_corner,rho,eta,reg_param]=l_curve(U,sm,Ct,'Tikh');
% Solve R(t) by Tikhonov regulrization:
[Rt,rho,eta] = tikhonov(U,sm,X,Ct,sqrt(reg_corner));

figure;plot(Rt);