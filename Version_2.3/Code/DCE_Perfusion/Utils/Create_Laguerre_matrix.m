function [ B ] = Create_Laguerre_matrix(knots,time_vec,alpha,poly_deg)
%Create_B_matrix Create the cubical B spline matrix

%   knots    - K ti's to set as knots
%   N        - number of total points
%   time_vec - time points for the function to be evaluated
%   poly_deg - polynomial degree
    
        
    K = max(size(knots));
    N = max(size(time_vec));
    M = poly_deg; % Poly degree
    
    B = zeros(N,K+M);

    for i = 1 : K+M %columns
            % Evaluates L{n, alpha}(x)
            B(:,i) = polyval(LaguerreGen(i, alpha), time_vec);
    end
    
end

  