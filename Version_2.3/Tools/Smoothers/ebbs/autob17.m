function [f,band]=autob17(x,y,h1a,h1b,p,order,msespan,M1,J1, ...
     J2,nterms,xgrid,varyfactor,ridgecoef,type,varest,coef,nskip,bandspan) ;
%
%  Copied from autob16.m                           
%                          
%
%  Last edited: 1/25/96
%
%  Automatic bandwidth for univariate local linear regression --- fits
%  a smooth curve to (x,y).
%
%  Usage: 
% [f,band]=autob17(x,y,h1a,h1b,p,order,msespan,M1,J1,J2,nterms,
%    xgrid,varyfactor,ridgecoef,type,varest,coef,nskip,bandspan);
%  
%  In calling program: x and y are the raw data
%                      h1a, h1b :  bandwidths between h1a and h1b will be tried
%                      
%                      order = order of derivative being estimated (= 0
%                              when the function itself is being estimated)
%                      p = degree of polynomials minus order ( = degree of
%                          polynomials when estimating the function itself)
%                      msespan = span for smoothing the mse---at the l-th
%                          point of xgrid, mse is averaged from the 
%                          (l-msespan)-th to (l+msespan)-th xgrid point 
%                          (except for truncation at boundaries)
%                      M1 = total number of fits for estimating bias at all
%                          bandwidths
%
%                      J1 = number of bandwidth points to the left
%                           used for fitting a curve to
%                           estimate bias at one bandwidth
%
%                      J2 = number of bandwidth points to the right
%                           used for fitting a curve to
%                           estimate bias at one bandwidth
%
%                      nterm = number of terms in the bias model
%                              The bias model is
%                                E(mhat(x;h)) = b_0 + b_1 h^(p+1) + ... +
%                                                b_nterms h^(p + nterms)
%
%                      xgrid     = m-point grid where f is estimated
%          
%                      coef = binheights/bincounts (use only if type = 1)
%
%                      type = type of smoothing problem
%
%                           = 1: estimation with a known variance function
%                             estimate, so the variance function
%                             estimate must be supplied as "varest"
%                             (type = 1 is required if order ~= 0)
%
%                           = 2: estimation with a known coefficient of
%                             variation of the response (e.g., 
%                             variance function estimation where the errors
%                             are in a scale family)
%                             ("coef" = (sigma/mean)^2 must be input)
%
%                           = 3: sigma^2/mean has a known value = "coef" 
%                             (e.g., in density estimation using binned data
%                             where the response
%                             is a known multiple of a Poisson variate)
%
%
%                      coef = coefficient of variation of the response 
%                             (when type = 2) or ratio of response variance
%                             to response mean (when type = 3).  Has no effect
%                             when type = 1
%
%                      varest = estimated variance function at each point 
%                               of xgrid (used only when
%                               type = 1)
%                     
%                      varyfactor = for binned data, the variance of the
%                        response after binning is varyfactor times the
%                        variance of the unbinned responses (this is a
%                        function of x) (for unbinned data use varyfactor = 1)
%
%                      ridgecoef = coefficient for local ridge regression
%                        (= 0 if local ridge regression is not wanted)
%
%                        varest = estimated variance function = estimated
%                         variance of the raw (unbinned) response as a 
%                         function of x (required in type = 1)
%
%                      nskip = number of grid points at boundaries that 
%                       are "skipped" when estimating MSE to get local
%                       bandwidths 
%
%                      bandspan = bandwidth for smoothing the bandwidth
%
%                
%  Returns:    f         = estimated function on xgrid (m by 1)
%              band      = bandwidth vector (gives bandwidth at each point on
%                          xgrid) (m by 1)
%
%
%      CALLS: lpolydb
%
%
%
%global msehat msevect msevect2 h2 h3 h4 fvect ;
h1a=log(h1a) ;
h1b=log(h1b) ;
m=size(xgrid,1) ;

if size(h1a,1) == 1;h1a=h1a*ones(m,1) ; end ;
if size(h1b,1) == 1;h1b=h1b*ones(m,1) ; end ;
h2 = zeros(m,M1) ;
for i=1:m ;

h2(i,:)=exp(linspace(h1a(i),h1b(i),M1) ) ;
end ;

	if size(varyfactor,1) == 1;
	varyfactor = varyfactor*ones(size(x,1),1) ;
	end ;
maxJ10 = max([0 J1]) ;
fvect = zeros(m,M1) ;
var2 = fvect ;
band = zeros(m,1) ;
msevect = zeros(m,M1-J2-maxJ10) ;
msehat = zeros(m,1) ;
   for k = 1:M1 ;
   [f,var]=lpolydb(x,y,h2(:,k),xgrid,order,p,varyfactor,ridgecoef) ;
   fvect(:,k) = f ;
   var2(:,k) = var ;
   end ;

H = ones(J2 + J1 + 1,nterms) ;
   for l = 1:m ;
       for ll = ( 1 + maxJ10 ):(M1-J2) ;
       bot = ll-J1 ;
       top = ll+J2 ;
         for j=2:(1+nterms)
         H(:,j) = h2(l,bot:top)'.^(p+j-1); %p = degree of poly - order of der.
         end ;
         GGTTT = H'*H;
         mmTT  = size(GGTTT,1);
       biascoef = inv(GGTTT + (.0000001 .* eye(mmTT)) ) * ( H'* fvect(l,(bot:top))' ) ;
       	   h3= h2(l,ll).^( (p+1):(p+nterms) ) ;     
       biassqu = ( biascoef(2:(nterms+1))'*h3' ).^2 ;  
           
	   if type == 1 ;	  
           msevect(l,ll-maxJ10) = biassqu + varest(l).*var2(l,ll) ;
	   elseif type == 2 ;
           msevect(l,ll-maxJ10) = biassqu + coef * (fvect(l,ll).^2) * var2(l,ll) ;
           elseif type == 3;
	   msevect(l,ll-maxJ10) = biassqu + coef*fvect(l,ll)*var2(l,ll) ;
	  	   end ;
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

[f,varfinal]=lpolydb(x,y,band,xgrid,order,p,varyfactor,ridgecoef) ;
  
