function [ Cb, Ce ] = Larsson_Differential_Equation( time_vec, Vb, Ve, Ca, Hct, F, Ki)
%Larsson_Differential_Equation - Drives the differential equation


idx_test = 2;
Cb   = zeros(size(time_vec));
Ce   = zeros(size(time_vec));

for tj = time_vec(2:end)
    
    % Initial condition
    if (idx_test==2)
        Cb(idx_test) = (1/Vb) * ( time_vec(idx_test) - time_vec(idx_test-1) ) * F * Ca(idx_test-1) ;
        Ce(idx_test) = 0;
        idx_test = idx_test + 1;
        continue;
    end
    
    % Regular case
    Cb(idx_test) = Cb(idx_test-1) + (1/Vb) * ( time_vec(idx_test) - time_vec(idx_test-1) )...
        * ( ( F * Ca(idx_test-1)   ) +...
        ( Ki*(1-Hct)* Ce(idx_test-1) ) -...
        ( Ki*Cb(idx_test-1) ) -...
        ( F*Cb(idx_test-1) ) ...
        );
    Ce(idx_test) = Ce(idx_test-1) + (1/Ve) * ( time_vec(idx_test) - time_vec(idx_test-1) )...
        * ( ( Ki* Cb(idx_test-1) ) -...
        ( Ki*(1-Hct)*Ce(idx_test-1) ) ...
        );
    
    idx_test = idx_test + 1;
    
end



end

