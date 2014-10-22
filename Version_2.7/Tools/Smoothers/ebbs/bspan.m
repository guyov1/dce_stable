function bsp = bspan(xgrid,x,p);
% 
%   Computes the bandwidth at each point of xgrid corresponding to
%    a span of p with x the vector of x values.
%
%  Last edited:  3/20/97
%
% USAGE: bsp = bspan(xgrid,x,p);
%
%       Copyright: David Ruppert
%
n = size(xgrid,1);
nx = size(x,1) ;
bsp = zeros(n,1) ;
np = ceil(nx*min([p;1])) ;


        for i = 1:n ;
        x2 = sort(abs(x-xgrid(i))) ;
        bsp(i) = x2(np)*max([p;1]) ;
        end ;

