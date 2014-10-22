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

% - finalK=5 means to cluster assuming mixture order of 5
% - omtr contains the cluster parameters of the mixture of order 5
disp('clustering data assuming optimal order of 5...');
[mtrs,omtr] = GaussianMixture(pixels, 20, 5, true);
for i=1:omtr.K
   disp(sprintf('\tCluster %d:', i));
   disp(sprintf('\t\tpi: %f', omtr.cluster(i).pb));
   disp([sprintf('\t\tmean: '), mat2str(omtr.cluster(i).mu',6)]);
   disp([sprintf('\t\tcovar: '), mat2str(omtr.cluster(i).R,6)]);
end
input('press <Enter> to continue...');


