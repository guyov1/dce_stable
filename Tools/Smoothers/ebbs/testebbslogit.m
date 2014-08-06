%	Program to test ebbslogit
%
%


n = 2000 ;

x= linspace(0,1,n)' ;
y = (rand(n,1) < .5 + .5*sin(10*x) ) ;

[f,xgrid,band]= ebbslogit(x,y, .2, .9,2) ;	%	p=2 instead of the
						%	default of 1

plot(xgrid,f) ;
set(gca,'Ylim',[0 1] );
