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
% load('/data/home/gilad/Desktop/Code/TmpRCTCEMat.mat');
load('/data/home/gilad/Desktop/Code/Tmp.mat');
% SimData=CurRCTCE(9,:);
SimData=TT;
N=numel(SimData);
Xs=1:N;
figure(9494);clf;
% plot(Xs,SimData,'*');
NSmoothers=4;
% ksr
gsubplot(NSmoothers,1);
r=ksr(Xs,SimData,20,N);
plot(Xs,SimData,'co',r.x,r.f,'r--','linewidth',2)
legend('data','regression','location','northwest');
title('Gaussian kernel regression')
% ksrlin
gsubplot(NSmoothers,2);
r=ksrlin(Xs,SimData,20,N);
plot(Xs,SimData,'co',r.x,r.f,'r--','linewidth',2)
legend('data','regression','location','northwest');
title('Gaussian kernel regression')
% lowess
[a b]=max(SimData);

gsubplot(NSmoothers,3);
   datain = [Xs;SimData]';
   f = 0.08;
   xdata = (1:0.2:N)';
   [dataout lowerLimit upperLimit xy] = lowess(datain,f,0,'',xdata);
   % A second plot to illustrate what values would be used, and how it
   % compares to the original function used to generate the data
plot(Xs,SimData,'co',xdata,xy(:,2),'r.','linewidth',2)
legend('data','regression','location','northwest');
title('lowess')
% supsmu
gsubplot(NSmoothers,4);
%   smo = supsmu(Xs,SimData,'Alpha',0);

Weights=SimData*0+1;
Weights(1:b)=linspace(1,100,b);
  smo = supsmu(Xs,SimData,'Alpha',0,'Weights',Weights);
plot(Xs,SimData,'co',Xs,smo,'r--','linewidth',2)
legend('data','regression','location','northwest');
title('supsmu')
