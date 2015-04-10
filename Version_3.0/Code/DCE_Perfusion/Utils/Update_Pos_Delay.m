function [ output_vec ] = Update_Pos_Delay( delay, input_vec )

%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

output_vec                     = zeros(size(input_vec));

output_vec(1         : delay ) = input_vec(end - delay + 1 : end        );
output_vec(delay + 1 : end   ) = input_vec(1               : end - delay);

end

