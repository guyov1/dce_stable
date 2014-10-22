function [phat,band]=ebbsdensity(x,ha,hb,p,order,span,M1,J1,J2,nterms,xgrid, ...
               ridgecoef) ;
%
%  Estimate a probability density with local polynomials and EBBS bandwidths
%
%  Last edited:  11/22/95
%
%
%     Estimates the density on xgrid
%
%     x = sample  (dim = n)
%
%     ha, hb = bandwidth range for smoothing 
%
%     p = degree of the local polynomials for y's - order
%
%     xgrid = locations where var. fn. is estimated
%
%                      order = order of derivative being estimated (= 0
%                              when the function itself is being estimated)
%
%                      span = span for smoothing the mse---at the l-th
%                          point of xgrid, mse is averaged from the 
%                          (l-msespan)-th to (l+msespan)-th xgrid point 
%                          (except for truncation at boundaries) and later
%                          for smoothing the bandwidth (span = 3 recommended)
%
%                      M1 = total number of fits for estimating bias at all
%                          bandwidths  (M1 = 12 is recommended)
%
%                      J1 = number of bandwidth points to the left
%                           used for fitting a curve to
%                           estimate bias at one bandwidth  (J1=1 recommended)
%
%                      J2 = number of bandwidth points to the right
%                           used for fitting a curve to
%                           estimate bias at one bandwidth (J2 = 1 recommended)
%
%                      nterm = number of terms in the bias model
%                              The bias model is
%                                E(mhat(x;h)) = b_0 + b_1 h^(p+1) + ... +
%                                                b_nterms h^(p + nterms)
%				(nterm = 1 (or 2) recommended)
%
%                      ridgecoef = ridge regression coefficient (= 0 for no
%                                  ridge regression)
%
%
%     This implementation uses the epan. kernel with support (-1,1).
%
%     returns:
%             phat = estimated order-th der of the density
%             band = bandwidths
%
%
%     calls the following functions:  lpolyres, autob17
%
%
if size(xgrid,1) == 1;
xgrid = xgrid' ;      % make xgrid a column
end ;
n = size(x,1) ;
nbin = 500 ;
minx=min(xgrid) ;
maxx=max(xgrid) ;
xgridsmall = linspace(minx,maxx,25)' ;  % small grid for estimating bandwidths
                                      % quickly
bincenters=linspace(minx,maxx,nbin)';
y = hist(x,bincenters)' ;             % bin data
coef = nbin/((maxx-minx)*n) ;
y = coef*y ;                          % rescale to be a prob. density
type = 3 ;                            % set autob17 parameters by default
coef = 0 ;
varyfactor = 1 ;
var= 1 ;
nskip = 1 ;
%   Get bandwidth on small grid
[phat,band]=autob17(bincenters,y,ha,hb,p,order,span,M1,J1, ...
     J2,nterms,xgridsmall,varyfactor,ridgecoef,type,var,coef,nskip,span) ;

band = interp1(xgridsmall,band,xgrid,'linear') ; % interpolate bandwidth to
                                                 % xgrid

[phat,var]= lpolydb(bincenters,y,band,xgrid, ...
            order,p,varyfactor,ridgecoef) ;     % get density on xgrid

phat = max([phat' ; zeros(1,size(phat,1))])' ;  % truncate at 0

phat = phat ./ ( (maxx-minx)*mean(phat) ) ;  % renormal to be a density
