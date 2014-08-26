function [h_t_est] = Wiener_Filter( x_t, y_t, Fs )
%Wiener_Filter Implements Wiener deconvolution
%   The model is y(t) = h(t)*x(t) + n(t)
%   The function gets y[num_voxels,nT] and dt*x[nT].
%   The result is estimation of h_t[nT].

% Set # of FFT points according to data size
num_fft_points = size(y_t,2);
num_voxels     = size(y_t,1);
sec_interval   = 1 / Fs;
min_interval   = sec_interval / 60;

% Frequency vector
f = Fs/2*linspace(-1,1,num_fft_points);

% Calculate FFT of x_y and y_t
Y_f = fft(y_t,num_fft_points,2);
X_f = fft(x_t,num_fft_points);


%% Rectangle estimation for SNR(f)

% Make the rectangle big enough comparing to |X(f)|^2
Max_val = max(X_f .* conj(X_f)) * 10^6 * 100 * 500000;
% The index will be chosen according to -0.05 to +0.05 Hz
% This way the bandwith will be 0.2 hertz which will capture up to
% 5 seconds change in signal
idx         = find(fftshift(f) > 0.1,1);

% Initiate SNR estimation
SNR_psd_est = zeros(size(X_f));

% Set the rectangle
% SNR_psd_est(1:idx)       = Max_val;
% SNR_psd_est(end-idx:end) = Max_val;

% Make it a full rectangle
SNR_psd_est(1:end) = Max_val;


%% Use Wiener Filter

% Create Wiener filter
G_f =  ( 1./ X_f ) .* ( ( X_f .* conj(X_f) ) ./ ( ( X_f .* conj(X_f) ) + ( 1./SNR_psd_est ) ) );

% Pass Y(f) through G(f) to estimate H(f)
H_f_est = repmat(G_f,num_voxels,1) .* Y_f;

% Convert H(f) to h(t) using IFFT
h_t_est = real(ifft(H_f_est,num_fft_points,2));

% Zero negative values
h_t_est(h_t_est<0) = 0;

end

