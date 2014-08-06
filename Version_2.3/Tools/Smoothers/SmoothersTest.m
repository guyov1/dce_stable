% cd('/u/peptibase3-ext/libermg1/Code');

P=path;
rmpath(genpath('/data/SPM5/'));
addpath(genpath('/data/spm8/'));
addpath('/data/Gilad/Base/Code');
addpath('/data/Gilad/SVN/');
addpath('/data/Gilad/SVN/DataInfrastructure');
addpath(genpath('/data/Gilad/Matlab stuff'))
addpath('/data/home/gilad/Desktop/T2/Code');

rmpath(genpath('/data/SVN'));
addpath(genpath('/data/Gilad/SVN'));

addpath(genpath('/data/home/gilad/Desktop/Code/Smoothers'))
%%
N=1000;
m=2;
b=3;
Sigma=50;
Xs=1:N;
SimData=Xs*m+b+randn(1,N)*Sigma;
figure(9494);clf;plot(Xs,SimData,'*');
%% ksr
r=ksr(Xs,SimData,20,N);
plot(Xs,SimData,'co',r.x,r.f,'r--','linewidth',2)
legend('data','regression','location','northwest');
title('Gaussian kernel regression')
%% ksrlin
r=ksrlin(Xs,SimData,20,N);
plot(Xs,SimData,'co',r.x,r.f,'r--','linewidth',2)
legend('data','regression','location','northwest');
title('Gaussian kernel regression')
%% lowess
   datain = [Xs;SimData]';
   f = 0.15;
   wantplot = 1;
   xdata = Xs';
   [dataout lowerLimit upperLimit xy] = lowess(datain,f,0);
   % A second plot to illustrate what values would be used, and how it
   % compares to the original function used to generate the data
plot(Xs,SimData,'co',Xs,dataout(:,3),'r--','linewidth',2)
legend('data','regression','location','northwest');
title('lowess')
%% supsmu
  smo = supsmu(Xs,SimData);
plot(Xs,SimData,'co',Xs,smo,'r--','linewidth',2)
legend('data','regression','location','northwest');
title('supsmu')
