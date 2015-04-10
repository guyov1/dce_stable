function [ X_mat ] = build_LQ_model_matrix( t, k, p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Initiate
X_mat = zeros(length(p), 3);

X_mat(:,1) = 1;

for i = 1 : p
    
   if i > k
       % Linear part
       X_mat(i,2) =   t(i) - t(k);
       % Quadratic part
       X_mat(i,3) = ( t(i) - t(k) ) ^2;
   else
       X_mat(i,2) = 0;
       % Quadratic part
       X_mat(i,3) = 0;
   end
   
    
end

end

