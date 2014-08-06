function [ est_gauss_filter_Wiener_noise,est_larss_filter_Wiener_noise, idx_fig  ] = Estimate_ht_Wiener( Sim_Struct, Verbosity, iter_num, avg_num, idx_fig)

% Take from struct variables used in local function
gauss_filter                  = Sim_Struct.gauss_filter;
larss_filter                  = Sim_Struct.larss_filter;
plot_flag                     = Sim_Struct.plot_flag;
Ignore_Gaussian_Calculation   = Sim_Struct.Ignore_Gaussian_Calculation;
min_interval                  = Sim_Struct.min_interval;
time_vec_minutes              = Sim_Struct.time_vec_minutes;
Sim_AIF_with_noise            = Sim_Struct.Sim_AIF_with_noise;
noise_to_add_gauss            = Sim_Struct.noise_to_add_gauss;
noise_to_add_larss            = Sim_Struct.noise_to_add_larss;
Sim_Ct_larss_kernel           = Sim_Struct.Sim_Ct_larss_kernel;
Sim_Ct_gauss_kernel           = Sim_Struct.Sim_Ct_gauss_kernel_noise;
Fs                            = Sim_Struct.Fs;

if strcmp(Verbosity,'Full')
    display('-I- Estimating h(t) using Wiener de-convolution...');
end

% Second, for Larsson filter
[f_larss, SNR_real_psd_larss, SNR_est_psd_larss, Y_f_larss, H_psd_larss, N_psd_larss,...
    Est_Wiener_Filter_real_PSD_larss, Est_Wiener_Filter_est_PSD_larss] = Wiener_Filter_4_Simulation( ...
    Sim_Ct_larss_kernel(:,iter_num,avg_num), noise_to_add_larss(:,iter_num,avg_num),...
    min_interval(iter_num)*Sim_AIF_with_noise(:,iter_num,avg_num), larss_filter(:,iter_num), Fs(iter_num));


% Tranpose the estimated filters by Wiener filter
est_larss_filter_Wiener_noise = transpose(Est_Wiener_Filter_est_PSD_larss);

% Once for gaussian filter
if (~Ignore_Gaussian_Calculation)
    
    [f_gauss, SNR_real_psd_gauss, SNR_est_psd_gauss, Y_f_gauss, H_psd_gauss, N_psd_gauss,...
        Est_Wiener_Filter_real_PSD_gauss, Est_Wiener_Filter_est_PSD_gauss] = Wiener_Filter_4_Simulation( ...
        Sim_Ct_gauss_kernel(:,iter_num,avg_num), noise_to_add_gauss(:,iter_num,avg_num),...
        min_interval(iter_num)*Sim_AIF_with_noise(:,iter_num,avg_num), gauss_filter(:,iter_num), Fs(iter_num));
    
    % Tranpose the estimated filters by Wiener filter
    est_gauss_filter_Wiener_noise = transpose(Est_Wiener_Filter_est_PSD_gauss);
else
    est_gauss_filter_Wiener_noise = zeros(size(est_larss_filter_Wiener_noise));
    
end


% Plotting frequency data
if (plot_flag)
    % PDF report title
    idx_fig    = idx_fig + 1;
    idx_string = ['idx_' num2str(idx_fig,'%03i')];
    AddToLog('Run_Output/',idx_string,'\\subsection*{\\underline{Wiener Filter}}');
    
    if (~Ignore_Gaussian_Calculation)
        
        filter_name = 'Gaussian';
        
        fig_num = figure;
        subplot(3,1,1);
        
        hold on;
        plot(f_gauss,abs(fftshift(Y_f_gauss))  ,'g*');
        plot(f_gauss,abs(fftshift(H_psd_gauss)),'yo');
        plot(f_gauss,abs(fftshift(N_psd_gauss)),'bd');
        hold off;
        
        title(['FFT Ct. PSD of H,N. Filter used: ' filter_name]);
        legend('Ct(f)','H psd(f)','N psd(f)');
        xlabel('Frequency (Hz)');
        
        subplot(3,1,2);
        
        hold on;
        plot(f_gauss,abs(fftshift(SNR_real_psd_gauss)),'g');
        plot(f_gauss,abs(fftshift(SNR_est_psd_gauss)),'r');
        hold off;
        title('SNR(f)');
        legend('SNR - real psd','SNR - est psd');
        xlabel('Frequency (Hz)');
        ylabel('|SNR(f)|')
        
        subplot(3,1,3);
        
        plot(time_vec_minutes,gauss_filter,'r*',time_vec_minutes,Est_Wiener_Filter_real_PSD_gauss,'gd',time_vec_minutes,Est_Wiener_Filter_est_PSD_gauss,'bo');
        title('Original vs. Estimated h(t)','FontWeight','bold');
        legend('Original h(t)','Estimated h(t)-PSD','Estimated h(t)-est PSD');
        xlabel('Time [Min]');
        
        % Print result to PDF
        file_name = ['Wiener_FFT_PSD_' filter_name '.png'];
        [idx_fig] = Print2Pdf(fig_num, idx_fig, file_name, './Run_Output/','', 'EstFilters');
        
        
    end
    
    
    filter_name = 'Larsson';
    
    fig_num = figure;
    subplot(3,1,1);
    
    hold on;
    plot(f_larss,abs(fftshift(Y_f_larss))  ,'g*');
    plot(f_larss,abs(fftshift(H_psd_larss)),'yo');
    plot(f_larss,abs(fftshift(N_psd_larss)),'bd');
    hold off;
    
    title(['FFT Ct. PSD of H,N. Filter used: ' filter_name]);
    legend('Ct(f)','H psd(f)','N psd(f)');
    xlabel('Frequency (Hz)');
    
    subplot(3,1,2);
    
    hold on;
    plot(f_larss,abs(fftshift(SNR_real_psd_larss)),'g');
    plot(f_larss,abs(fftshift(SNR_est_psd_larss)),'r');
    hold off;
    title('SNR(f)');
    legend('SNR - real psd','SNR - est psd');
    xlabel('Frequency (Hz)');
    ylabel('|SNR(f)|')
    
    subplot(3,1,3);
    
    plot(time_vec_minutes,larss_filter,'r*',time_vec_minutes,Est_Wiener_Filter_real_PSD_larss,'gd',time_vec_minutes,Est_Wiener_Filter_est_PSD_larss,'bo');
    title('Original vs. Estimated h(t)','FontWeight','bold');
    legend('Original h(t)','Estimated h(t)-PSD','Estimated h(t)-est PSD');
    xlabel('Time [Min]');
    
    % Print result to PDF
    file_name = ['Wiener_FFT_PSD_' filter_name '.png'];
    [idx_fig] = Print2Pdf(fig_num, idx_fig, file_name, './Run_Output/','', 'EstFilters');
    
end

end

