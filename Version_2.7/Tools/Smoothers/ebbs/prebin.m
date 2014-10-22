function [x2,y2,varyfactor,varest] = prebin(x,y,nbin) ;
%
%
%      Last edited 6/16/95
%
%  Bins (x,y) data using local linear regression with "boxcar" weights
%
%  The (x,y) data are sorted according to their x-values and put into
%  nonoverlapping groups of approximately "nx" data points.  For the ith
%  group 
%    x2(i) = mean of the x's
%    y2(i) = estimated mean of the regression function = intercept at x2(i)
%    varest(i) = MSE from this fit
%    varyfactor = factor to multiply by varest(i) to get the estimated
%                 variance of y2(i)
%   
%
%   Also returns cvvarest = estimated global coefficient of variation of the
%                           varest values
%
%   INPUT: x,y = raw data
%          nbin = number of bins
%
%   The program computes n=dim(x) and nx= floor(n/nbin)
%
%   USAGE: function [x2,y2,varyfactor,varest] = prebin(x,y,nbin) ;
% 
%
global nx ;
[x,I]=sort(x) ;
y=y(I) ;
n= rows(x) ;
nx = floor(n/nbin) ; 

x2=zeros(nbin,1) ;
y2 = x2 ;
varyfactor = x2 ;
varest = x2 ;

        for m = 0:(nbin-1) ;
            if m < nbin-1 ;
            x3 = [ ones(nx,1) x((m*nx+1):((m+1)*nx)) ] ;
            y3 = y((m*nx+1):((m+1)*nx)) ;
            else;
            x3=[ ones(n-m*nx,1) x((m*nx+1:n)) ] ;
            y3 = y((m*nx+1):n) ;
            end ;
        xmean = mean(x3(:,2)) ;
        x2(m+1) = xmean ;
        x3(:,2) = x3(:,2) - xmean ;
                ixx = inv(x3'*x3) ;
        beta = ixx*(x3'*y3) ;
        y2(m+1) = beta(1) ;
        varest(m+1) =sum((y3-x3*beta).^2)/(rows(x3)-2) ;
        varyfactor(m+1) = ixx(1,1) ;
        end ;


