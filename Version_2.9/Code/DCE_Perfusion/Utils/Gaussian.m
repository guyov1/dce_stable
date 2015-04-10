
function to_fit = Gaussian(t,mean,var,amplitude)
    to_fit = amplitude * (1/sqrt(2*pi*(var))).*exp(- ( (t - mean).^2 / (2*var) ) ); 
end
