function [Cost, reversedParams] = ReverseTransformation(Func,transformedParams,LB,UB, time_vec, estF)

reversedParams = zeros(size(transformedParams));

for i = 1 : length(reversedParams)
    bias                  = ( LB(i) + UB(i) ) / 2 ;
    scale                 =   UB(i) - bias;
    
    %transformed_params(i) = ( atan( params(i) ) / ( pi/2 ) )*scale   + bias; % param is bounded to low_limit-upper_limit
    reversedParams(i)     = tan( ((transformedParams(i) - bias ) / scale ) * (pi/2) );
    
end

est_F_noise = estF;
Cost        = Func(reversedParams, time_vec);

end