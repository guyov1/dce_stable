function [m,bandmean,var,bandvar]=ebbsmvparam(x,y,hm,hmeana,hmeanb, ...
     hvara,hvarb,pmean,pvar,order,xgrid,msespan,M1,J1,J2,nterms, ...
     varyfactor,ridgecoef,nskip,bandspan) ;
%
%  Estimate mean and variance with EBBS bandwidths
%
%  Last edited:  3/12/98
%
%
%     Estimates the local mean and variance on xgrid
%
%     Calls these subroutines: ebbsvar, autob17, lpolyvar, lpolyres, lpolydb,
%	prebin
%
if size(xgrid,1) == 1;
xgrid = xgrid' ;
end ;

[var,bandvar]=ebbsvar(x,y,hm,hvara,hvarb,pvar, ...
         xgrid,msespan,M1,J1,J2,nterms,varyfactor,ridgecoef,nskip,bandspan) ;

type = 1 ;
coef = 0 ;
[m,bandmean]=autob17(x,y,hmeana,hmeanb,pmean,order,msespan,M1,J1, ...
     J2,nterms,xgrid,varyfactor,ridgecoef,type,var,coef,nskip,bandspan) ;
