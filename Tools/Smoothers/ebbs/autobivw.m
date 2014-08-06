function [f,band,xgrid1,xgrid2]=autobivw(x,y,h1a,h1b,p,order,msespan,M1,J1, ...
     J2,nterms,gridsize1,gridsize2, ...
     varyfactor,ridgecoef,type,varest,coef,bandspan,weights) ;
%
%  Copied from autobivw: allows weighting in addition to kernel weights 
%                         
%                          
%
%  Last edited: 1/11/95
%
%  Automatic bandwidth for bivariate local linear regression --- fits
%  a smooth curve to (x,y).
%
%  Usage:  [f,band,xgrid1,xgrid2]=autobiv(x,y,h1a,h1b,p,order,msespan,M1,J1, ...
%                   J2,nterms,gridsize1,gridsize2, ...
%                   varyfactor,ridgecoef,type,varest,coef,bandspan,weights) ;
%
%
%  
%  In calling program: x and y are the raw data
%                      h1a, h1b :  bandwidths between h1a and h1b will be tried
%                      
%                      order = order of derivative being estimated (= 0
%                              when the function itself is being estimated)
%                      p = degree of polynomials minus order ( = degree of
%                          polynomials when estimating the function itself)
%
%                      msespan = span for smoothing the mse
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
%      CALLS: lpolydb, rows
%
%
%
global msehat msevect msevect2 h2 h4 fvect min1 min2 max1 max2 xgrid h2 ...
 band2 band3 xgridb  xgrid1old xgrid2old xgrid1 xgrid2 ;

min1=min(x(:,1)) ;
min2=min(x(:,2)) ;
max1=max(x(:,1)) ;
max2=max(x(:,2)) ;

xgrid1= linspace(min1,max1,gridsize1)' ;
xgrid2= linspace(min2,max2,gridsize1)' ;

xgrid=[kron( ones(gridsize1,1), xgrid1) kron( xgrid2,ones(gridsize1,1)) ] ;
% this appears to be an error - detected 12/7

xgrid=[kron(xgrid1,ones(gridsize1,1)) kron(ones(gridsize1,1),xgrid2) ] ; 
% added 12/7

m=rows(xgrid) ;
h1a=log(h1a) ;
h1b=log(h1b) ;
h2=exp(linspace(h1a,h1b,M1) ) ;
	if rows(varyfactor) == 1;
	varyfactor = varyfactor*ones(rows(x),1) ;
	end ;
maxJ10 = max([0 J1]) ;
fvect = zeros(m,M1) ;
var2 = fvect ;
band = zeros(m,1) ;

msevect = zeros(m,M1-J2-maxJ10) ;
msehat = zeros(m,1) ;
   for k = 1:M1 ;

   [f,var]=lpolybivw(x,y,h2(:,k),xgrid,order,p,varyfactor,ridgecoef,0,weights) ;
   fvect(:,k) = f ;
   var2(:,k) = var ;
   end ;

H = ones(J2 + J1 + 1,nterms) ;
   for l = 1:m ;
       for ll = ( 1 + maxJ10 ):(M1-J2) ;
       bot = ll-J1 ;
       top = ll+J2 ;
         for j=2:(1+nterms)
         H(:,j) = h2(bot:top)'.^(p+j-1); %p = degree of poly - order of der.
         end ;
       biascoef = inv(H'*H) * ( H'* fvect(l,(bot:top))' ) ;
       	   h4= h2(ll).^( (p+1):(p+nterms) ) ;     
       biassqu = ( biascoef(2:(nterms+1))'*h4' ).^2 ;  
           
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
h4 =linspace(h2(1+maxJ10)+2*eps,h2(M1-J2)-2*eps,n2) ;
	
   for l=1:m ;
   msevect2(l,:) = interp1(h2((maxJ10+1):(M1-J2)),msevect(l,:),h4,'cubic')' ;
   end ;

kopt = 2 ;



   for l = 1:m ;
   weight = epan( (xgrid(:,1) - xgrid(l,1))./msespan ) ...
            .*epan( (xgrid(:,2) - xgrid(l,2))./msespan ) ;
   weight = weight ./ sum(weight) ;
   mse = weight' * msevect2 ;
  
   im=  n2;
	for lll=1:(n2-kopt) ;
          if min(  mse( (lll+1:lll+kopt) )  ) > mse(lll) ;
               if im == n2 ;
               im =lll;
               end ;
	  end;
	end ;
   band(l) = h4(im) ; 
   msehat(l) = mse(im) ;
   end ;



	if bandspan > 0 ;   %  SMOOTH THE BANDWIDTH
        band2=band ;
        	for l=1:m ;
		weight = epan( (xgrid(:,1) - xgrid(l,1))./bandspan ) ...
		    .* epan( (xgrid(:,2) - xgrid(l,2))./bandspan ) ;
		weight = weight ./ sum(weight) ;
		band(l)= weight' * band2 ;
		end ;
	band = band2 ;
	end;


xgrid1old=xgrid1;
xgrid2old=xgrid2;

xgrid1 = linspace(min1,max1,gridsize2)' ;
xgrid2 = linspace(min2,max2,gridsize2)' ;

xgrid1b=kron( xgrid1, ones(gridsize2,1) ) ;
xgrid2b=kron( ones(gridsize2,1), xgrid2);



xgridb =[xgrid1b xgrid2b] ;

band2 = reshape(band,gridsize1,gridsize1) ;

band2 = interp2(xgrid1old',xgrid2old,band2,xgrid1',xgrid2) ;

band = reshape(band2,gridsize2^2,1) ;

[f,var]=lpolybivw(x,y,band,xgridb,order,p,varyfactor,ridgecoef,0,weights) ;

