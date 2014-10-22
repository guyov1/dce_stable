
load('AIFmat.mat');
load('Ct.mat');

lambda=1;
% We want first derivative to be negative (R is a descending function of
% time)
% D=zeros(size(AIFmat));
zero_vec=zeros(size(Ct));

D_main_diagonal=(-1)*ones(size(Ct));
D_second_diagonal=ones(length(Ct)-1,1);

D1=diag(D_main_diagonal,0);
D2=diag(D_second_diagonal,1);

D_temp=D1+D2;

D=(D_temp+lambda*(D_temp*D_temp));

Rt = lsqlin(AIFmat,Ct,D,zero_vec);

err=norm(AIFmat*Rt-Ct);

figure;plot(Rt);title(['error = ',num2str(err)],'FontSize',12);

figure;plot(abs((D*D)*Rt)-mean(abs((D*D)*Rt)));
