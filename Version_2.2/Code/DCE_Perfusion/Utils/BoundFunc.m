function Cost=BoundFunc(Func,param1,param2,LB,UB)

x=(atan(param1)/(pi/2))+5; % X is bounded to 4-6
y=(atan(param1)/(pi/2))*3+1; % Y is bounded -2 to 4

Cost=Func(x,y);


% ----------------------
% 
% OtherSpaceBestParams=fminsearch(BoundFunc);
% BestParams=ReverseTransformation(OtherSpaceBestParams);






% ToDo:
% Optimizations - always try
%     1. lsqcurvefit (trust-region-reflective)
%     2. lsqcurvefit (LM)
%     3. fminsearch with cost RMS (?, or better is you have idea, like weighted)
%     4. patternsearch
%     
% Finer search: After all the stuff, let F as free variable and see if it improves
% 
% Do a graph checking the cost as function of F, to see if we can find a more
% correct F starting from the estimated one, or if we need to estimate anyway


% lsqcurvefit(fun,x0,xdata,ydata)
% 
% RMSCost=@(vec) sum(vec.^2);
% 
% CostFunc=@(x) RMSCost(fun(x,xdata)-ydata)
% 
% fminsearch(CostFunc,x0)