% Double Gaussian function
function to_fit = DoubleGaussian(t,mean_1,var_1,amplitude_1,mean_2,var_2,amplitude_2)
    to_fit = amplitude_1 * (1/sqrt(2*pi*(var_1))).*exp(- ( (t - mean_1).^2 / (2*var_1) ) ) + ...
             amplitude_2 * (1/sqrt(2*pi*(var_2))).*exp(- ( (t - mean_2).^2 / (2*var_2) ) ); 
end
