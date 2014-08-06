function fit = ebbsmvdef(x,y,order,msespan,M1,J1,J2,bandspan,nterms,...
	pmean,pvar,varyfactor,ridgecoef) ;
%
%    Last edited:  6/10/2000
%
%
%  Calls ebbsmvparam with default values of
%  USAGE:  fit = ebbsmvdef(x,y,order,msespan,M1,J1,J2,nterms,...
	varyfactor,ridgecoef,nskip,pmean,pvar) ;
%
%
%		INPUT (required)
%  	x and y are the independend and dependent variable (n by 1)
%
%		INPUT (optional)
%  	order = order of derivative being estimated (default = 0)
%
%
%		OUTPUT (all in a structure named "fit")
%  	mean = estimated derivative of mean function
%  	bandmean = bandwidths used for mean
%  	var = estimated variance function
%  	bandvar = bandwidth for var
%  	xgridmean = grid that mean is calculated on
%  	xgridvar = grid that var is calculated on
%
%  Calls these subroutines:  autob17, lpolydb, ebbsvar, lpolyres, bspan
%                            ebbsmvparam, lpolyvar, prebin
%
%	Copyright: David Ruppert

if exist('param.order') == 0 ;
	order = 0 ;
else ;
	order = param.order ;
end ;

if exist('param.msespan') == 0 ;
	msespan = 0 ;
else ;
	msespan = param.msespan ;
end ;


if exist('param.M1') == 0 ;
	M1 = 14 ;
else ;
	M1 = param.M1 ;
end ;

if exist('param.J1') == 0 ;
	J1 = 1 ;
else ;
	J1 = param.J1 ;
end ;

if exist('param.J2') == 0 ;
	J2 = 2 ;
else ;
	J2 = param.J2 ;
end ;

if exist('param.bandspan') == 0 ;
	bandspan =4 ;
else ;
	bandspan = param.bandspan ;
end ;

if exist('param.nterms') == 0 ;
	nterms = 2 ;
else ;
	nterms = param.nterms ;
end ;

if exist('param.pmean') == 0 ;
	pmean = 2 ;
else ;
	pmean = param.pmean ;
end ;

if exist('param.pvar') == 0 ;
	pvar = 0 ;
else ;
	pvar = param.pvar ;
end ;


if exist('param.xxx') == 0 ;
	xxx = 0 ;
else ;
	xxx = param.xxx ;
end ;


if exist('param.varyfactor') == 0 ;
	varyfactor = 0 ;
else ;
	varyfactor = param.varyfactor ;
end ;

if exist('param.ridgecoef') == 0 ;
	ridgecoef = 0 ;
else ;
	ridgecoef = param.ridgecoef ;
end ;


minx=min(x) ;
maxx = max(x) ;

n = size(y,1) ; 
xgrid=linspace(minx,maxx,25)';
hm = bspan(x,x,8/n) ;
hmeana = bspan(xgrid,x,10/n) ;
hmeanb = bspan(xgrid,x,1.2) ;
hvara = bspan(xgrid,x,20/n) ;
hvarb = bspan(xgrid,x,1.5) ;
eval(string) ;


 [mean,bandmean,var,bandvar]=ebbsmvparam(x,y,hm,hmeana,hmeanb, ...
     hvara,hvarb,pmean,pvar,order,xgrid,msespan,M1,J1,J2,nterms,varyfactor, ...
     ridgecoef,nskip,bandspan) ;

xgridmean = linspace(minx,maxx,200)';
bandmean = interp1(xgrid,bandmean,xgridmean,'linear') ;

[mean,var2] = lpolydb(x,y,bandmean,xgridmean, ...
            order,pmean,varyfactor,ridgecoef) ;     % get mean on xgrid2


xgridvar = xgrid ;

fit = struct('mean',mean,'xgridmean',xgridmean,'bandmean', ...
	bandmean,'var',var,'xgridvar',xgridvar,'bandvar',bandvar) ;
