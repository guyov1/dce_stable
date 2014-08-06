function fit = ebbsmv01_raymodified(x,y,param) ;
%
%    	Last edited:  6/12/2000
%
%	Copied from ebbsmvdef
%
%  	Calls ebbsmvparam with default values of tuning parameters,
%		though non-default values can be specified in
%		argument "param" which is a structure
%  	USAGE:  fit = ebbsmvdef(x,y,param) ;
%
%	Example: fit = ebbsmv01(x,y,struct('minx',0,'maxx',1)) ;
%		In this example, all parameters are at their
%		default values except "minx" and "maxx".
%
%		INPUT (required)
%  	x and y are the independent and dependent variable (n by 1)
%	param
%		INPUT (optional - put as fields into the structure
%				called "param")
%	
%  	order = order of derivative being estimated (default = 0)
%	msespan = span for smoothing mse (default = 0 )
%	M1 (default = 14)
%	J1 (default = 1)
%	J2 (default = 2)
%	bandspan (default = 4)
%	nterms (default = 2)
%	pmean (default = 2)
%	pvar (default = 0)
%	varyfactor (default = 1)
%	ridgecoef (default = 0)
%	nskip (default = 0)
%	minx = left endpoint of xgridmean and xgridvar (default = min(x))
%	maxx = right endpoint of xgridmean and xgridvar (default = max(x))
%
%		OUTPUT (all in a structure named "fit")
%  	mean = estimated derivative of mean function
%  	bandmean = bandwidths used for mean
%  	var = estimated variance function
%  	bandvar = bandwidth for var
%  	xgridmean = grid that mean is calculated on
%  	xgridvar = grid that var is calculated on
%
%  	Calls these subroutines:  autob17, lpolydb, ebbsvar, lpolyres, bspan
%		ebbsmvparam, lpolyvar, prebin
%
%	Copyright: David Ruppert

if nargin < 3 ;
	param = '' ;
end ;

if isfield(param,'order') == 0 ;
	order = 0 ;
else ;
	order = param.order ;
end ;

if isfield(param,'msespan') == 0 ;
	msespan = 0 ;
else ;
	msespan = param.msespan ;
end ;


if isfield(param,'M1') == 0 ;
	M1 = 14 ;
else ;
	M1 = param.M1 ;
end ;

if isfield(param,'J1') == 0 ;
	J1 = 1 ;
else ;
	J1 = param.J1 ;
end ;

if isfield(param,'J2') == 0 ;
	J2 = 2 ;
else ;
	J2 = param.J2 ;
end ;

if isfield(param,'bandspan') == 0 ;
	bandspan =4 ;
else ;
	bandspan = param.bandspan ;
end ;

if isfield(param,'nterms') == 0 ;
	nterms = 2 ;
else ;
	nterms = param.nterms ;
end ;

if isfield(param,'pmean') == 0 ;
	pmean = 2 ;
else ;
	pmean = param.pmean ;
end ;

if isfield(param,'pvar') == 0 ;
	pvar = 0 ;
else ;
	pvar = param.pvar ;
end ;



if isfield(param,'varyfactor') == 0 ;
	varyfactor = 1 ;
else ;
	varyfactor = param.varyfactor ;
end ;

if isfield(param,'ridgecoef') == 0 ;
	ridgecoef = 0 ;
else ;
	ridgecoef = param.ridgecoef ;
end ;

if isfield(param,'nskip') == 0 ;
	nskip = 0 ;
else ;
	nskip = param.nskip ;
end ;

if isfield(param,'minx') == 0 ;
	minx=min(x) ;
else ;
	minx = param.minx ;
end ;

if isfield(param,'nxgrid') == 0 ;
	nxgrid=101 ;
else ;
	nxgrid = param.nxgrid ;
end ;

if isfield(param,'maxx') == 0 ;
	maxx = max(x) ;
else ;
	maxx = param.maxx ;
end ;


n = size(y,1) ; 
xgrid=linspace(minx,maxx,25)';
hm = bspan(x,x,6/n) ;
hmeana = bspan(xgrid,x,10/n) ;
hmeanb = bspan(xgrid,x,1.2) ;
hvara = bspan(xgrid,x,20/n) ;
hvarb = bspan(xgrid,x,1.5) ;



% The idea is to first fit the function on a small
% grid, and then expand to the extire vector
% It's at the latter point that I have required that
% the X's first be sorted into ascending order
[mean,bandmean,var,bandvar]=ebbsmvparam(x,y,hm,hmeana,hmeanb, ...
     hvara,hvarb,pmean,pvar,order,xgrid,msespan,M1,J1,J2,nterms,varyfactor, ...
     ridgecoef,nskip,bandspan) ;

 xgridmean = linspace(minx,maxx,nxgrid)';
bandmean = interp1(xgrid,bandmean,xgridmean,'linear') ;
[mean,var2] = lpolydb(x,y,bandmean,xgridmean, ...
            order,pmean,varyfactor,ridgecoef) ;     % get mean on xgrid2

all_xgridmean = x;
all_bandmean = interp1(xgridmean,bandmean,x,'linear') ;
all_mean     = interp1(xgridmean,mean,x,'linear');
%[all_mean,all_var2] = lpolydb(x,y,all_bandmean,all_xgridmean, ...
%            order,pmean,varyfactor,ridgecoef) ;     % get mean on xgrid2


xgridvar = xgrid ;
fit = struct('mean',mean,'xgridmean',xgridmean,'bandmean', ...
	bandmean,'var',var,'xgridvar',xgridvar,'bandvar',bandvar, ...
	     'all_mean',all_mean) ;
