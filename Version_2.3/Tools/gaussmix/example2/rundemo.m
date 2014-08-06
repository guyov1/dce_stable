clear all;

if exist('GaussianMixture')~=2
   pathtool;
   error('the directory containing the Cluster program must be added to the search path');
end

disp('generating data...');
mkdata;
clear all;
traindata1 = load('TrainingData1');
traindata2 = load('TrainingData2');
testdata = load('TestingData');
input('press <Enter> to continue...');
disp(' ');

% [mtrs, omtr] = GaussianMixture(pixels, initK, finalK, verbose)
% - pixels is a NxM matrix, containing N training vectors, each of M-dimensional
% - start with initK=20 initial clusters
% - finalK=0 means estimate the optimal order
% - verbose=true displays clustering information
% - mtrs is an array of structures, each containing the cluster parameters of the
%   mixture of a particular order
% - omtr is a structure containing the cluster parameters of the mixture with
%   the estimated optimal order
disp('clustering class 1...');
[mtrs,class1] = GaussianMixture(traindata1, 20, 0, false);
disp(sprintf('\toptimal order K*: %d', class1.K));
for i=1:class1.K
   disp(sprintf('\tCluster %d:', i));
   disp(sprintf('\t\tpi: %f', class1.cluster(i).pb));
   disp([sprintf('\t\tmean: '), mat2str(class1.cluster(i).mu',6)]);
   disp([sprintf('\t\tcovar: '), mat2str(class1.cluster(i).R,6)]);
end
input('press <Enter> to continue...');
disp(' ');

disp('clustering class 2...');
[mtrs,class2] = GaussianMixture(traindata2, 20, 0, false);
disp(sprintf('\toptimal order K*: %d', class2.K));
for i=1:class2.K
   disp(sprintf('\tCluster %d:', i));
   disp(sprintf('\t\tpi: %f', class2.cluster(i).pb));
   disp([sprintf('\t\tmean: '), mat2str(class2.cluster(i).mu',6)]);
   disp([sprintf('\t\tcovar: '), mat2str(class2.cluster(i).R,6)]);
end
input('press <Enter> to continue...');
disp(' ');

disp('performing maximum likelihood classification...');
disp('for each test vector, the following calculates the log-likelihood given each of the two classes, and classify');
disp('the first half of the samples are generated from class 1, the remaining half from class 2');
disp(' ');
likelihood=zeros(size(testdata,1), 2);
likelihood(:,1) = GMClassLikelihood(class1, testdata);
likelihood(:,2) = GMClassLikelihood(class2, testdata);
class=ones(size(testdata,1),1);
class(find(likelihood(:,1)<=likelihood(:,2)))=2;
for n=1:size(testdata,1)
   disp([mat2str(testdata(n,:),4), sprintf('\tlikelihood: '), mat2str(likelihood(n,:),4), sprintf('\tclass: %d',  class(n))]);
end




