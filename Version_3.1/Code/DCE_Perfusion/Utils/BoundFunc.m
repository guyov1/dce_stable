function [boundedParams] = BoundFunc(unboundedParams,LB,UB)

boundedParams = zeros(size(unboundedParams));

for i = 1 : length(unboundedParams)
    bias                  = ( LB(i) + UB(i) ) / 2 ;
    scale                 =   UB(i) - bias;
    boundedParams(i) = ( atan( unboundedParams(i) ) / ( pi/2 ) )*scale   + bias; % param is bounded to low_limit-upper_limit
end

end


% function [Cost, transformed_params] = BoundFunc(Func,params,LB,UB, time_vec, estF)
% 
% transformed_params = zeros(size(params));
% 
% for i = 1 : length(params)
%     bias                  = ( LB(i) + UB(i) ) / 2 ;
%     scale                 =   UB(i) - bias;
%     transformed_params(i) = ( atan( params(i) ) / ( pi/2 ) )*scale   + bias; % param is bounded to low_limit-upper_limit
% end
% 
% est_F_noise = estF;
% Cost        = Func(transformed_params, time_vec);
% 
% end

% ----------------------
% 
% OtherSpaceBestParams = fminsearch(BoundFunc);
% BestParams           = ReverseTransformation(OtherSpaceBestParams);


% Usage exp.
% Func   = @(x) sum(x.^2)
% RMSCost  = @(vec) sum(vec.^2); 
% CostFunc = @(x) RMSCost( fun(x,xdata) - ydata)

% LB       = [0 0 0];
% UB       = [5 10 5];
% params   = [0.5 7 3];

% [Cost, transformedParams] = BoundFunc(CostFunc, params, LB, UB);

% best_transformedParams    = fminsearch(CostFunc,x0)
% reversedParams            = ReverseTransformation(transformedParams, LB, UB)


% ToDo:
% Optimizations - always try
%     1. lsqcurvefit (trust-region-reflective)
%     2. lsqcurvefit (LM)
%     3. fminsearch with cost RMS (?, or better if you have idea, like weighted)
%     4. patternsearch
%     
% Finer search: After all the stuff, let F be a free variable and see if it improves
% 
% Do a graph checking the cost as function of F, to see if we can find a more
% correct F starting from the estimated one, or if we need to estimate anyway


% lsqcurvefit(fun,x0,xdata,ydata)
% 
% RMSCost  =@(vec) sum(vec.^2); 
% CostFunc =@(x) RMSCost(fun(x,xdata)-ydata)
% 
% fminsearch(CostFunc,x0)