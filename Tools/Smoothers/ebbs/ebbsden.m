function [phat,band,xgrid]=ebbsden(x,xgrid,p,ha,hb,order) ;
%
%  Estimate a probability density with local polynomials and EBBS bandwidths
%
%  Uses default parameters in ebbsdensity
%
%  Last edited:  12/16/97
%
%
%     Estimates the density on xgrid
%
%		REQUIRED INPUT
%    	 x = sample  (dim = n)
%
%		OPTIONAL INPUT
%     	xgrid = locations where density is estimated (default 
%			linspace(min(x)-.1*range,max(x)+.1,range,80))
%
%     	ha = lower limit of bandwidth range (default = range/20)
%
%    	hb = upper limit of bandwidth range (default = .8*range)
%
%     	p = degree of the local polynomials for y's - order (default = 1)
%
%     	order = order of derivative being estimated (= 0
%            when the function itself is being estimated) (default = 0)
%
%
%     		OUTPUT
%	phat = estimated order-th der of the density
% 	band = bandwidths
%	xgrid = same as input if xgrid is input.  Otherwise, same as default.
%
%     This implementation uses the epan. kernel with support (-1,1).
%
% Call is [phat,band,xgrid]=ebbsden(x,xgrid,p,ha,hb,order);
%
%     calls the following functions:  lpolydb, autob17, 
%                                     ebbsdensity
%
%
%
%
rangex = range(x) ;

if nargin < 2 ;
xgrid = linspace(min(x)-.1*rangex,max(x)+.1*rangex,100) ;
end ;

if nargin < 3 ;
p = 1 ;
end ;

if nargin < 6 ;
order = 0 ;
end ;


if nargin < 5 ;
	if p == 1 ;
	ha = rangex/20 ;
	hb = .8*rangex ;
	elseif p == 2 ;
	ha = rangex/15 ; hb = rangex ;
	elseif p > 2 ;
	ha = rangex/10 ;hb = 1.2*rangex ;
	end ;
end ;


span = 5 ;
M1 = 12 ;
J1 = 1 ;
J2 = 2 ;
nterms = 2 ;
ridgecoef = 0;
[phat,band]=ebbsdensity(x,ha,hb,p,order,span,M1,J1,J2,nterms,xgrid, ...
               ridgecoef) ;

