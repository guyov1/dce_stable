function [f,xgrid,band]=ebbslogit(x,y,h1a,h1b,p,msespan,M1,J1, ...
     J2,nterms,xgrid,nskip,bandspan) ;
%
%  Copied from autob17.m                           
%                          
%
%  Last edited: 9/7/98
%  Currently setup only for local logistic regression
%
%  EBBS bandwidth for univariate local logistic regression --- fits
%  a smooth curve to (x,y).
%
%  Usage: 
% [f,xgrid,band]=ebbslogit(x,y,h1a,h1b,p,msespan,M1,J1,J2,nterms,
%    xgrid,nskip,bandspan);
%
%  				INPUT (required)
%	x and y are the raw data
%	h1a, h1b :  bandwidths between h1a and h1b will be tried
%                      
%				INPUT (optional)
%	p = degree of polynomials (dafault = 1)
%
%	msespan = span for smoothing the mse---at the l-th
%	point of xgrid, mse is averaged from the 
%	(l-msespan)-th to (l+msespan)-th xgrid point 
%	(except for truncation at boundaries) (default = 0)
%
%	M1 = total number of fits for estimating bias at all
%	bandwidths (default = 14)
%
%	J1 = number of bandwidth points to the left
%	used for fitting a curve to estimate bias at one bandwidth
%	(default = 1)
%
%	J2 = number of bandwidth points to the right
%	used for fitting a curve to estimate bias at one bandwidth
%	(default = 2)
%
%	nterm = number of terms in the bias model
%	The bias model is E(mhat(x;h)) = b_0 + b_1 h^(p+1) + ... +
%	b_nterms h^(p + nterms) (default = 2)
%
%	xgrid is initially a m-point grid where bandwidth is estimated.
%	At the end, xgrid is a 5m-point grid where the function
%	is evaluated.  (The EBBS bandwidth is interpolated to the
%	finer bandwidth.)
%	(default for initial xgrid = linspace(min(x),max(x),20)')
%
%	nskip = number of grid points at boundaries that 
%	are "skipped" when estimating MSE to get local
%       bandwidths (default = 0) 
%
%	bandspan = bandwidth for smoothing the bandwidth
%	(default = 4)
%
%                
%  Returns:    f         = estimated function on xgrid (m by 1)
%              band      = bandwidth vector (gives bandwidth at each point on
%                          xgrid) (m by 1)
%
%
%      CALLS: locallogit
%
%
%
if nargin < 13 ;
	bandspan = 4 ;
end ;

if nargin < 12 ;
	nskip = 0 ;
end ;

if nargin < 11 ;
	xgrid = linspace(min(x),max(x),20)' ;
end ;

if nargin < 10 ;
	nterms = 2 ;
end ;

if nargin < 9 ;
	J2=2 ;
end ;

if nargin < 8 ;
	J1 =1 ;
end ;

if nargin < 7 ;
	M1 = 14 ;
end ;

if nargin < 6 ;
	msespan = 0 ;
end ;

if nargin < 5 ;
	p = 1 ;
end ;


h1a=log(h1a) ;
h1b=log(h1b) ;
m=size(xgrid,1) ;

if size(h1a,1) == 1;h1a=h1a*ones(m,1) ; end ;
if size(h1b,1) == 1;h1b=h1b*ones(m,1) ; end ;
for i=1:m ;
h2(i,:)=exp(linspace(h1a(i),h1b(i),M1) ) ;
end ;
	;
maxJ10 = max([0 J1]) ;
fvect = zeros(m,M1) ;
var2 = fvect ;
band = zeros(m,1) ;
msevect = zeros(m,M1-J2-maxJ10) ;
msehat = zeros(m,1) ;
   for k = 1:M1 ;
   [f,var]=locallogit(x,y,h2(:,k),xgrid,p) ;
   fvect(:,k) = f ;
   var2(:,k) = var ;
   end ;

H = ones(J2 + J1 + 1,nterms) ;
   for l = 1:m ;
       for ll = ( 1 + maxJ10 ):(M1-J2) ;
       bot = ll-J1 ;
       top = ll+J2 ;
         for j=2:(1+nterms)
         H(:,j) = h2(l,bot:top)'.^(p+j-1); % p = degree of poly
         end ;
       biascoef = inv(H'*H) * ( H'* fvect(l,(bot:top))' ) ;
       	   h3= h2(l,ll).^( (p+1):(p+nterms) ) ;     
       biassqu = ( biascoef(2:(nterms+1))'*h3' ).^2 ;  
       msevect(l,ll-maxJ10) = biassqu + var(l) ;
       end ;
   end ;

n2=35 ;
msevect2=ones(m,n2) ;              %  INTERPOLATE MSE TO FINER GRID
	
h4=zeros(m,n2) ;

   for l=1:m ;
   h4(l,:) =linspace(h2(l,1+maxJ10)+2*eps,h2(l,M1-J2)-2*eps,n2) ;
  
msevect2(l,:) = interp1(h2(l,(maxJ10+1):(M1-J2)),msevect(l,:),h4(l,:), ...
                   'cubic') ;

   end ;

kopt = 2 ;
mseweight = [(1:msespan+1) fliplr(1:msespan) ] ;  %  WEIGHTS FOR SMOOTHING
                                                  %  THE MSE
mseweight = mseweight ./ sum(mseweight) ; 

   for l = (1+msespan):(m-msespan) ;
   index = (l-msespan):(l+msespan) ;
   mse = mseweight*msevect2(index,:) ;
   im=  n2;
	for lll=1:(n2-kopt) ;
          if min(  mse( (lll+1:lll+kopt) )  ) > mse(lll) ;
               if im == n2 ;
               im =lll;
               end ;
	  end;
	end ;
   band(l) = h4(l,im) ; 
   msehat(l) = mse(im) ;
   end ;

msns = msespan + nskip ;  %  ESTIMATE MSE NEAR BOUNDARIES
	if msns > 0 ;
	band(1:msns)   = band(msns +1)*ones(msns,1) ;
	band((m-msns+1):m) = band(m-msns)*ones(msns,1) ;
	end ;


	if bandspan > 0 ;   %  SMOOTH THE BANDWIDTH
        band2=band ;
	delta = bandspan.*(xgrid(2) - xgrid(1)) ;
		for i = 1:m;
		wt = max(delta - abs(xgrid-xgrid(i)),zeros(m,1)) ;
		wt = wt ./ sum(wt) ;
		band2(i) = wt'*band ;
		end  ;
	band = band2 ;
	end;
xgridold = xgrid ;
xgrid = linspace(min(xgrid),max(xgrid),5*m)' ;
band = interp1(xgridold,band,xgrid,'linear') ;
[f,var]=locallogit(x,y,band,xgrid,p) ;
  
