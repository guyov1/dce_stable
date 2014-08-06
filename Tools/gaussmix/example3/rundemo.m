clear all;

if exist('GaussianMixture')~=2
   pathtool;
   error('the directory containing the Cluster program must be added to the search path');
end

disp('generating data...');
mkdata;
clear all;
pixels = load('data');
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
disp('estimating optimal order and clustering data...');
[mtrs,omtr] = GaussianMixture(pixels, 20, 0, true);
disp(sprintf('\toptimal order K*: %d', omtr.K));
for i=1:omtr.K
   disp(sprintf('\tCluster %d:', i));
   disp(sprintf('\t\tpi: %f', omtr.cluster(i).pb));
   disp([sprintf('\t\tmean: '), mat2str(omtr.cluster(i).mu',6)]);
   disp([sprintf('\t\tcovar: '), mat2str(omtr.cluster(i).R,6)]);
end
input('press <Enter> to continue...');
disp(' ');

% with omtr containing the optimal clustering with order K*, the following
% split omtr into K* classes, each containing one of the subclusters of
% omtr
disp('split the optimal clustering into classes each containing a subcluster...');
mtrs = SplitClasses(omtr);
input('press <Enter> to continue...');
disp(' ');

disp('performing maximum likelihood classification...');
disp('for each test vector, the following calculates the log-likelihood given each of the classes, and classify');
disp(' ');
likelihood=zeros(size(pixels,1), length(mtrs));
for k=1:length(mtrs)
   likelihood(:,k) = GMClassLikelihood(mtrs(k), pixels);
end
class=ones(size(pixels,1),1);
for k=1:length(mtrs)
   class(find(likelihood(:,k)==max(likelihood,[],2)))=k;
end
for n=1:size(pixels,1)
   str = [mat2str(pixels(n,:),4), sprintf('\tlikelihood: '), mat2str(likelihood(n,:),4), sprintf('\tclass: %d',  class(n))];
   disp(str);
end
