% Wiener filter for simulation
function [f, SNR_real_psd, SNR_est_psd, Y_f, H_psd, N_psd, h_t_est_real_PSD, h_t_est_est_PSD] = Wiener_Filter_4_Simulation( y_t, n_t, x_t, h_t, Fs )

%Wiener_Deconvolution Implements Wiener deconvolution
%   The model is y(t) = h(t)*x(t) + n(t)
%   The function gets y[nT], n[nT] and dt*x[nT].
%   The result is estimation of h_t[nT].

% Number of FFT points will be the size of data vector
num_fft_points = max(size(y_t,1));

% Prepare frequency vector
f = Fs/2*linspace(-1,1,num_fft_points);

% Calculate FFT of relevant vectors
Y_f = fft(y_t,num_fft_points);
X_f = fft(x_t,num_fft_points);
N_f = fft(n_t,num_fft_points);
H_f = fft(h_t,num_fft_points);

% Calculate H's and N's PSD for SNR calculation
Hxx    = abs(H_f).^2/length(h_t)/Fs;
Nxx    = abs(N_f).^2/length(n_t)/Fs;

%H_psd  = dspdata.psd(Hxx,'SpectrumType','Twosided','Fs',Fs);
%H_psd  = H_psd.data;
H_psd = pwelch(Hxx, hanning(length(h_t)),0, length(h_t), Fs, 'twosided');
%figure; plot(H_psd,'or'); hold on; plot(H_psd_2,'*g'); hold off;

%N_psd  = dspdata.psd(Nxx,'SpectrumType','Twosided','Fs',Fs);
%N_psd  = N_psd.data;
N_psd = pwelch(Nxx, hanning(length(n_t)),0, length(n_t), Fs, 'twosided');
%figure; plot(N_psd,'or'); hold on; plot(N_psd_2,'*g'); hold off;


% Calculating SNR
SNR_real_psd = transpose(H_psd) ./ transpose(N_psd);

% Estimate psd SNR by place where max value drops by 100
SNR_est_psd = zeros(size(SNR_real_psd));
%idx         = find( abs(SNR_real_psd) < max(abs(SNR_real_psd)) / 100,1);

% Rectangle - make it big enough comparing to |X(f)|^2
Max_val = max(X_f .* conj(X_f)) * 10^6 * 100 * 500000;

% The index will be chosen according to -0.1 to +0.1 Hz
% This way the bandwith will be 0.2 hertz which will capture up to 5 seconds change in signal
idx                      = find(fftshift(f) > 0.1,1);
SNR_est_psd(1:idx)       = Max_val;
SNR_est_psd(end-idx:end) = Max_val;

% If time resolution is better than 1 seconds interval, 
% assume SNR(f) of f > 1Hz is noise
% if (Fs > 1)
%     idx                      = find(fftshift(f) > 0.5,1);
%     SNR_est_psd(1:idx)       = Max_val;
%     SNR_est_psd(end-idx:end) = Max_val;
% else
%     SNR_est_psd(1:end)       = Max_val;
% end


% Create Wiener filter
G_f_real_PSD =  ( 1./ X_f ) .* ( ( X_f .* conj(X_f) ) ./ ( ( X_f .* conj(X_f) ) + ( 1./SNR_real_psd' ) ) );
G_f_est_PSD  =  ( 1./ X_f ) .* ( ( X_f .* conj(X_f) ) ./ ( ( X_f .* conj(X_f) ) + ( 1./SNR_est_psd' ) ) );

% Pass Y(f) through G(f) to estimate H(f)
H_f_est_real_PSD = G_f_real_PSD .* Y_f;
H_f_est_est_PSD  = G_f_est_PSD  .* Y_f;

% Convert H(f) to h(t) using IFFT
h_t_est_real_PSD = real(ifft(H_f_est_real_PSD,num_fft_points));
h_t_est_est_PSD  = real(ifft(H_f_est_est_PSD,num_fft_points));

% Zero negative values
h_t_est_real_PSD(h_t_est_real_PSD<0) = 0;
h_t_est_est_PSD(h_t_est_est_PSD<0) = 0;

end

