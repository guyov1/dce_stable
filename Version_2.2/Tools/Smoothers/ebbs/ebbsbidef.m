function [mhat,band,varmhat,xgrid1,xgrid2]= ebbsbidef(x,y,covy,string) ;

% Bivariate EBBS with default parameters


p=2 ;         %  Start setting defaults


p1 = .1 ;  % For spans
p2 = 1 ;
order = 0 ;
msespan = 0.001 ;
ridgecoef = 0 ;
nskip = 0 ;
bandspan = 4 ;     %  Finished setting defaults
xgrid = 0 ;
eval(string) ;     %  Change some of the defaults as specified by "string"

str = '';

order = 1 ;

[mhat,band,varmhat,xgrid1,xgrid2]=autobic(x,y,xgrid,p1,p2,p,order,msespan, ...
             covy,bandspan,str) ;
