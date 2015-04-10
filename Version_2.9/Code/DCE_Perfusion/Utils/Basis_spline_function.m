function [ output ] = Basis_spline_function(x, i, m, tau, M)
%Basis_spline_function Calculate Basis Spline functions
%   Input -
%           x   - the point in which to calculate the basis function
%           i   - knot number (j in first equation in wiki)
%           m   - degree of polynomial (k in wiki)
%           tau - knots list (tau_1, tau_2, ...)
%
%    Implemented according to the appendix named "Computational Considerations for Splines"
%    in the book "The Elements of Statistical Learning".


% Pad the knots with 2*m-1 zeros
%padding_zeros = zeros(1,2*m);
%tau_padded = [tau padding_zeros];

% Pad the end with the last value (2*m times)
padding_knots_end  = repmat(tau(end),1,M);
padding_knots_init = repmat(tau(1)  ,1,M);
tau_padded = [padding_knots_init tau padding_knots_end];

%padding_knots  = repmat(tau(end),1,2*m);
%tau_padded = [tau padding_knots];


if (m == 1) % Stop point
    
%     % Convention to avoid dividing by zero
%     if ( tau_padded(i) == tau_padded(i+1) )
%         output = 0;
%         return;
%     end
    
    %
    %     if ( i+1 > length(tau_padded) )
    %         output = 1;
    %         return;
    %     end
    
    if ( (x >= tau_padded(i)) && (x < tau_padded(i+1)) )
        output = 1;
        return;
        
    else
        output = 0;
        return;
    end
    
else % Recursion
    
    %     % Convention to avoid dividing by zero
    %     if ( all( tau_padded(i:i+m) == tau_padded(i) ) )
    %        output = 0;
    %        return;
    %     elseif ( ( (tau_padded(i+m)-tau_padded(i+1)) == 0 ) || ((tau_padded(i+m-1)-tau_padded(i)) == 0) )
    %        output = 0;
    %        return;
    %     end
    %
    %     output = ( ( (x-tau_padded(i))   / (tau_padded(i+m-1)-tau_padded(i)) )*Basis_spline_function(x,i,m-1,tau,M)   )+...
    %              ( ( (tau_padded(i+m)-x) / (tau_padded(i+m)-tau_padded(i+1)) )*Basis_spline_function(x,i+1,m-1,tau,M) ) ;
    %     return;
    
    
    %     % Convention to avoid dividing by zero
    %     if ( tau_padded(i+m) == tau_padded(i)  )
    %        output = 0;
    %        return;
    %     end
    
    
% %     if ( i+m > length(tau_padded) )
% %         alpha_i_m   = (x-tau_padded(i))   / (tau_padded(end)-tau_padded(i));
% %     else
% %         alpha_i_m   = (x-tau_padded(i))   / (tau_padded(i+m)-tau_padded(i));
% %     end
% %     
% %     
% %     if ( i+1+m > length(tau_padded) )
% %         if (i+1>length(tau_padded))
% %             alpha_i_1_m = (x-tau_padded(end))   / (tau_padded(end)-tau_padded(end));
% %         else
% %             alpha_i_1_m = (x-tau_padded(i+1))   / (tau_padded(end)-tau_padded(i+1));
% %         end
% %         
% %     else
% %         alpha_i_1_m = (x-tau_padded(i+1))   / (tau_padded(i+1+m)-tau_padded(i+1));
% %         
% %     end
% % 
% %     
% %     % Convention to avoid dividing by zero
% %     if ( i+m > length(tau_padded) ||tau_padded(i+m) == tau_padded(i)  )
% %         alpha_i_m = 0;
% %     elseif ( (x==tau_padded(i))  && (tau_padded(i+m)==tau_padded(i)) )
% %         alpha_i_m=0;
% %         
% %     end
% %     
% %     % Convention to avoid dividing by zero
% %     if ( ( i+1+m > length(tau_padded) ) || ( tau_padded(i+1+m) == tau_padded(i+1) ) )
% %         alpha_i_1_m = 0;
% %     elseif (  (x==tau_padded(i+1)) && (tau_padded(i+1+m)==tau_padded(i+1)) )
% %         alpha_i_1_m = 0;
% %     end
% %     
% %     
% %     output = alpha_i_m*Basis_spline_function(x,i,m-1,tau,M) +...
% %         ( 1-alpha_i_1_m )*Basis_spline_function(x,i+1,m-1,tau,M)  ;
% %     return;
    

    part1 = 0;
    part2 = 0;
    if (tau_padded(i+m-1)-tau_padded(i) ~= 0 )
        part1 = (x-tau_padded(i))   / (tau_padded(i+m-1)-tau_padded(i))*Basis_spline_function(x,i,m-1,tau,M);
    end
    if (tau_padded(i+m)-tau_padded(i+1) ~= 0)
        part2 = (tau_padded(i+m) - x)   / (tau_padded(i+m)-tau_padded(i+1))*Basis_spline_function(x,i+1,m-1,tau,M);
    end
    
    output = part1 + part2;
    
%     output = (x-tau_padded(i))   / (tau_padded(i+m-1)-tau_padded(i))*Basis_spline_function(x,i,m-1,tau,M) +...
%              (tau_padded(i+m) - x)   / (tau_padded(i+m)-tau_padded(i+1))*Basis_spline_function(x,i+1,m-1,tau,M)  ;
    return;
    
    
end


end

% DEBUG

% close all;
% y1 = [];
% y2 = [];
% time_vec = 0:0.01:1;
% for x = time_vec
%
%     y1 = [y1 Basis_spline_function(x,1,3,[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0])];
%     y2 = [y2 Basis_spline_function(x,2,3,[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0])];
%
% end
% figure;hold on; plot(time_vec,y1,'*g');plot(time_vec,y2,'-r');hold off;