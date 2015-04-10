function [ output_vec ] = ReScale_AIF( scale_factor, input_vec )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

output_vec         = scale_factor * (input_vec / max(input_vec) ); %[mM]

end

