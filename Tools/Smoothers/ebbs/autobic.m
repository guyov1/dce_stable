function [mhat,band,varmhat,xgrid1,xgrid2]=autobic(x,y,xgrid,p1,p2,p,...
   order,msespan,covy,bandspan,str) ;
% 
%                         
%                          
%
%  Last edited: 10/11/96
%
%  Automatic bandwidth for bivariate local linear regression.  Uses lpolybic.m
%  to compute the local polynomial fits.  The bandwidths are selected on
%  an automatically generated grid and then interpolated onto "xgrid",
%  which can be user specified or automatically generated.  The final
%  estimate is computed on xgrid.
%
%   USAGE:
%   [mhat,band,varmhat,xgrid1,xgrid2]=autobic(x,y,xgrid,p1,p2,p,...
%   order,msespan,covy,bandspan) ;
%
%
%
%  
%  In calling program: x and y are the raw data
%
%                      p1, p2 :  spans between p1 and p2 will be tried
%                      
%     order = 1 estimate the function itself
%             2 estimate 1st der. wrt x1
%             3 estimate 1st der. wrt x2
%             4 estimate 2nd der. wrt x1
%             5 estimate 2nd der. wrt x1
%             4 estimate mixed 2nd der. wrt x1 and x2
%
%                      p = degree of polynomials 
%
%                      msespan = span for smoothing the mse
%
%                      M1 = total number of fits for estimating bias at all
%                        bandwidths (Default is 12 but can be changed by "str".)
%
%                      J1 = number of bandwidth points to the left
%                           used for fitting a curve to
%                           estimate bias at one bandwidth (Default is 1
%                           but can be changed by "str".)
%
%                      J2 = number of bandwidth points to the right
%                           used for fitting a curve to
%                           estimate bias at one bandwidth (Default is 2
%                           but can be changed by "str".)
%
%                      nterm = number of terms in the bias model
%                              The bias model is
%                                E(mhat(x;h)) = b_0 + b_1 h^(p+1) + ... +
%                                                b_nterms h^(p + nterms)
%                              (Default is 2 but can be changed by "str".)
%
%                      xgrid = m-point grid where f is estimated. If
%                         xgrid==0 is specified by the user, 
%                          then a rectangular grid of size
%                         (gridsize2)^2 is automatically generated and used
%                         as xgrid.  
%
%                        The bandwidths are first selected
%                        on a rectangular grid of size (gridsize1)^2 and then
%                        interpolated onto xgrid
%
%                         gridsize1 = 5 and gridsize2=20
%                         are defaults, but these can be changed by "str".
%                             
%
%                      ridgecoef = coefficient for local ridge regression
%                        (= 0 if local ridge regression is not wanted)
%                         The default is 0, but this can be changed by "str".
%
%
%
%                      bandspan = bandwidth for smoothing the bandwidth
%
%                
%  Returns:    mhat         = estimated function on xgrid (m by 1).  If
%                             xgrid is automatically generated, then mhat
%                             is (gridsize2)^2.
%
%              band      = bandwidth vector (gives bandwidth at each point on
%                          xgrid) (m by 1)
%
%              xgrid1 and xgrid2: If xgrid is automatically generated then
%              it is xgrid1 by xgrid2 (each of these is a gridsize2 vector)
%
%      If xgrid == 0 is specified, so that xgrid is automatically generated,
%         then one can plot mhat as follows:
%
%             surf(xgrid1,xgrid2,mhat) ; % Surface lighting plot
%  or
%             contour(xgrid1,xgrid2,mhat) ;
%
%
%
%      CALLS: lpolybic, rows, cols, epan
%
%
%

%  Set defaults
M1 = 12 ;
J1 = 1 ;
J2 = 2 ;
nterms = 2 ;
gridsize1 = 5 ;
gridsize2 = 20 ;
ridgecoef = 0 ;

% Update defaults
eval(str) ;
xgriduser = xgrid ;

min1=min(x(:,1)) ;  %  Create the grid for bandwidth selection
min2=min(x(:,2)) ;
max1=max(x(:,1)) ;
max2=max(x(:,2)) ;

xgrid1= linspace(min1,max1,gridsize1)' ;
xgrid2= linspace(min2,max2,gridsize1)' ;

xgrid=[kron(xgrid1,ones(gridsize1,1)) kron(ones(gridsize1,1),xgrid2) ] ;


h1a = bspanbi(xgrid,x,p1) ;
h1b = bspanbi(xgrid,x,p2) ;

m=rows(xgrid) ;
h1a=log(h1a) ;
h1b=log(h1b) ;

h2 = zeros(m,M1) ;
for i=1:m ;

h2(i,:)=exp(linspace(h1a(i),h1b(i),M1) ) ;
end ;

maxJ10 = max([0 J1]) ;

% Set up storage matrices

mhatvect = zeros(m,M1) ;  % For storage of mhat for various bandwidths
varmhatvect = mhatvect ;
band = zeros(m,1) ;
msevect = zeros(m,M1-J2-maxJ10) ;
msehat = zeros(m,1) ;

   for k = 1:M1 ;

   [mhat,varmhat]=lpolybic(x,y,h2(:,k),xgrid,order,p,covy,ridgecoef,0) ;
   mhatvect(:,k) = mhat ;
   varmhatvect(:,k) = varmhat ;
   end ;

H = ones(J2 + J1 + 1,nterms) ;
   for l = 1:m ;
       for ll = ( 1 + maxJ10 ):(M1-J2) ;
       bot = ll-J1 ;
       top = ll+J2 ;
         for j=2:(1+nterms)
         H(:,j) = h2(bot:top)'.^(p+j-1); %p = degree of poly - order of der.
         end ;
       biascoef = inv(H'*H) * ( H'* mhatvect(l,(bot:top))' ) ;
       	   h4= h2(ll).^( (p+1):(p+nterms) ) ;     
       biassqu = ( biascoef(2:(nterms+1))'*h4' ).^2 ;  
           
	   	  
           msevect(l,ll-maxJ10) = biassqu + varmhatvect(l,ll) ;
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
band2 = reshape(band,gridsize1,gridsize1) ;

	if xgriduser == 0 ;
	xgrid1 = linspace(min1,max1,gridsize2)'  ;
	xgrid2 = linspace(min2,max2,gridsize2)'  ;
	xgrid1b = kron( xgrid1, ones(gridsize2,1) ) ;
	xgrid2b = kron( ones(gridsize2,1), xgrid2) ;
	xgridb =[xgrid1b xgrid2b] ;

	else ;
	xgridb = xgriduser ;
	end ;

band =interp2(xgrid1old',xgrid2old,band2,xgridb(:,1),xgridb(:,2) ) ;
[mhat,varmhat] = lpolybic(x,y,band,xgridb,order,p,covy,ridgecoef,0) ;

if xgriduser == 0 ;
mhat = reshape(mhat,gridsize2,gridsize2) ;
end ;

