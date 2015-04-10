function [ x ] = Regularized_Sol( A, b, x0, lambda, norm_flag, pos_flag )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% Ridge Regression

% Normalize flag for ridge() function
 warning('off','MATLAB:rankDeficientMatrix');
if (norm_flag == 1)
    x = ridge(b, A,lambda,0);
    x = x(2:end);
else
    x = ridge(b,A,lambda,1);
end
warning('on','MATLAB:rankDeficientMatrix');

% Remove zeros
if pos_flag
    x(x<0) = 0;
end

%% Tichonov
%[U, s, V] = csvd(A,'full');
%x         = tikhonov(U, s, V, b, lambda, x0);

%[U, s, V] = csvd(A);
%x         = tikhonov(U, s, V, b, lambda);


end

