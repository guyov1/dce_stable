function [ B ] = Create_B_matrix(knots,time_vec,poly_deg)
%Create_B_matrix Create the cubical B spline matrix

%   knots    - K ti's to set as knots
%   N        - number of total points
%   time_vec - time points for the function to be evaluated
%   poly_deg - polynomial degree
    
    display('-I- Creating splines basis matrix...');
        
    K = max(size(knots));
    N = max(size(time_vec));
    M = poly_deg; % Poly degree
    
    B = zeros(N,K+M);

    for i = 1 : K+M %columns
        
        for j = 1 : N % rows
            B(j,i) = Basis_spline_function( time_vec(j), i, poly_deg, knots, poly_deg);
        end
    end
    
end

