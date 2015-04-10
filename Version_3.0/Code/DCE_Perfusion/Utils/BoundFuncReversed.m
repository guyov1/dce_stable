function [unboundeddParams] = BoundFuncReversed(boundedParams,LB,UB)

unboundeddParams = zeros(size(boundedParams));

for i = 1 : length(unboundeddParams)
    bias                  = ( LB(i) + UB(i) ) / 2 ;
    scale                 =   UB(i) - bias;
    
    %transformed_params(i) = ( atan( params(i) ) / ( pi/2 ) )*scale   + bias; % param is bounded to low_limit-upper_limit
    unboundeddParams(i)     = tan( ((boundedParams(i) - bias ) / scale ) * (pi/2) );
    
end


end