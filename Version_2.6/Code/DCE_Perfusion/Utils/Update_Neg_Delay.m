function [ output_vec ] = Update_Neg_Delay( delay, input_vec )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% Change to a positive index
delay = -1 * delay;

output_vec = zeros(size(input_vec));

output_vec(end - delay + 1 : end        ) = input_vec(1         : delay);
output_vec(1               : end - delay) = input_vec(delay + 1 : end);


end

