function [varsmooth,band]=ebbsvar(x,y,hm,h1a,h1b,pvar, ...
         xgrid,msespan,M1,J1,J2,nterms,varyfactor,ridgecoef,nskip,bandspan) ;

%
%
%     ebbsvar    --- For smoothing standardized squared residuals to estimate
%                    a variance function (or one of its derivatives)
%                    with automatic selection of the
%                    bandwidth
%
%     last edit:	June 13, 2000
%
%     Estimates the local variance (est = varest) on xgrid
%
%     x = independent variable  (dim = n)
%     y = dependent variable  (dim = n)
%     hm = bandwidth for estimating mean to get residuals
%     h1a, h1b = bandwidth range for smoothing squared residuals
%     pvar = degree of the local polynomials for smoothing squared residuals
%     xgrid = locations where var. fn. is estimated
%     nskip = number of xgrid points at boundaries that are "skipped"
%             when calculating MSE to get local bandwidths (nskip = 0 or1
%             is recommended)
%
%     This implementation uses the epan. kernel with support (-1,1).
%
%     returns:
%             varest = estimated order-th der of variance function
%             band = bandwidths for smoothing squared residuals
%
%
%     calls the following functions:  lpolyres, 
%                                     autob17, lpolyvar	
%
%     USAGE: function [varest,band]=ebbsvar(x,y,hm,h1a,h1b,p, ...
%         xgrid,msespan,M1,J1,J2,nterms,varyfactor,ridgecoef,nskip,bandspan) ;
%
%
global res2 resvar;
[res2,resvar]=lpolyres(x,y,hm,2,varyfactor,0) ; 	
					% 	Get residuals with
					%	an essentially bias-free 
					%	quadratic fit
stdres2 = res2 ./ resvar ;
nbin = floor(size(x,1)/10) ;
[x2,smstdres2,varyfactor2,varvarest] = prebin(x,stdres2,nbin) ;
coef = mean( varvarest ./ (smstdres2.*smstdres2) ) ;
type = 2 ;
varest2 = 0 ;
order = 0 ;


[varest,band] = autob17(x,stdres2,h1a,h1b,pvar,order,msespan, ...
	M1,J1,J2,nterms, ...
   	xgrid,varyfactor,ridgecoef,type,varest2,coef,nskip,bandspan) ;
varsmooth = lpolyvar(x,y,hm,band,xgrid,pvar,varyfactor,ridgecoef) ;

