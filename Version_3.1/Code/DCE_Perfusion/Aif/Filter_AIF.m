function [ Sim_Struct ] = Filter_AIF( Sim_Struct, Verbosity )

tic;

if ~strcmp(Verbosity,'None')
    display('-I- Starting filtering AIFs...');
end

% Take from struct variables used in local function
min_interval                     = Sim_Struct.min_interval;
time_vec_minutes                 = Sim_Struct.time_vec_minutes;
F                                = Sim_Struct.F;
Ve_larss                         = Sim_Struct.Ve_larss;
Vb_larss                         = Sim_Struct.Vb_larss;
Vd                               = Sim_Struct.Vd;
E                                = Sim_Struct.E;
Sim_AIF_delayed_no_noise         = Sim_Struct.Sim_AIF_delayed_no_noise;
Sim_AIF_HighRes_delayed_no_noise = Sim_Struct.Sim_AIF_HighRes_delayed_no_noise;
num_iterations                   = Sim_Struct.num_iterations;
num_averages                     = Sim_Struct.num_averages;
gauss_filter                     = Sim_Struct.gauss_filter;
Sim_AIF_no_noise                 = Sim_Struct.Sim_AIF_no_noise;
Sim_AIF_with_noise               = Sim_Struct.Sim_AIF_with_noise;
larss_filter                     = Sim_Struct.larss_filter;
idx_fig                          = Sim_Struct.idx_fig;
SNR_ratio                        = Sim_Struct.SNR_ratio;
plot_flag                        = Sim_Struct.plot_flag;
Ignore_Gaussian_Calculation      = Sim_Struct.Ignore_Gaussian_Calculation;
High_res_min                     = Sim_Struct.High_res_min;
gauss_filter_HighRes             = Sim_Struct.gauss_filter_HighRes;
larss_filter_HighRes             = Sim_Struct.larss_filter_HighRes;

% "filter" works for vectors only, hence, loop for each iteration

% Initiate matrices
Sim_Ct_gauss_kernel          = zeros(size(Sim_AIF_delayed_no_noise));
Sim_Ct_larss_kernel          = zeros(size(Sim_AIF_delayed_no_noise));
noise_to_add_gauss           = zeros(size(Sim_AIF_delayed_no_noise));
noise_to_add_larss           = zeros(size(Sim_AIF_delayed_no_noise));
Sim_Ct_larss_kernel_high_res = zeros(size(Sim_AIF_HighRes_delayed_no_noise));
noise_to_add_larss_high_res  = zeros(size(Sim_AIF_HighRes_delayed_no_noise));

% Filter each Ca(t) and add noise
[Sim_Ct_gauss_kernel, Sim_Ct_larss_kernel, Sim_Ct_larss_kernel_high_res, noise_to_add_gauss, noise_to_add_larss, noise_to_add_larss_high_res] = Filter_AIF_Process(Sim_Struct, Sim_Ct_gauss_kernel, Sim_Ct_larss_kernel, Sim_Ct_larss_kernel_high_res, noise_to_add_gauss, noise_to_add_larss, noise_to_add_larss_high_res);

Sim_Ct_gauss_kernel_noise          = Sim_Ct_gauss_kernel          + noise_to_add_gauss;
Sim_Ct_larss_kernel_noise          = Sim_Ct_larss_kernel          + noise_to_add_larss;
Sim_Ct_larss_kernel_noise_high_res = Sim_Ct_larss_kernel_high_res + noise_to_add_larss_high_res;

% Zero negative values (non realistic)
Sim_Ct_gauss_kernel_noise(Sim_Ct_gauss_kernel_noise<0) = 0;
Sim_Ct_larss_kernel_noise(Sim_Ct_larss_kernel_noise<0) = 0;

% Put result in struct
Sim_Struct.noise_to_add_gauss                 = noise_to_add_gauss;
Sim_Struct.noise_to_add_larss                 = noise_to_add_larss;
Sim_Struct.Sim_Ct_gauss_kernel_noise          = Sim_Ct_gauss_kernel_noise;
Sim_Struct.Sim_Ct_larss_kernel_noise          = Sim_Ct_larss_kernel_noise;
Sim_Struct.Sim_Ct_larss_kernel_noise_high_res = Sim_Ct_larss_kernel_noise_high_res;
Sim_Struct.Sim_Ct_gauss_kernel                = Sim_Ct_gauss_kernel;
Sim_Struct.Sim_Ct_larss_kernel                = Sim_Ct_larss_kernel;
Sim_Struct.Sim_Ct_larss_kernel_high_res       = Sim_Ct_larss_kernel_high_res;

% Iteration to display
iter_display = 1;

% Draw simulated AIF, filter, result and result + noise
if (plot_flag)
    
    fig_num = figure;
    subplot(3,2,1);
    hold on;
    %     h1 = plot(time_vec_minutes,Sim_AIF_no_noise(:,iter_display,iter_display),'r');
    %     h2 = plot(time_vec_minutes,Sim_AIF_no_noise(:,iter_display,iter_display),'r*');
    %     h3 = plot(time_vec_minutes,Sim_AIF_with_noise(:,iter_display,iter_display),'b');
    %     h4 = plot(time_vec_minutes,Sim_AIF_with_noise(:,iter_display,iter_display),'bd');
    h1 = plot(time_vec_minutes,Sim_AIF_delayed_no_noise(:,iter_display,iter_display),'r');
    h2 = plot(time_vec_minutes,Sim_AIF_delayed_no_noise(:,iter_display,iter_display),'r*');
    h3 = plot(time_vec_minutes,Sim_AIF_with_noise(:,iter_display,iter_display),'b');
    h4 = plot(time_vec_minutes,Sim_AIF_with_noise(:,iter_display,iter_display),'bd');
    title('Simulated AIF (possibly delayed)','FontWeight','bold');
    xlabel('Time [Min]');
    ylabel('[mM]');
    hold off;
    legend([h2 h4],'Sim AIF','Sim AIF+noise');
    
    if ~Ignore_Gaussian_Calculation
        subplot(3,2,3);
        plot(time_vec_minutes,gauss_filter(:,iter_display),time_vec_minutes,gauss_filter(:,iter_display),'*');
        title('Gaussian filter used','FontWeight','bold');
        xlabel('Time [Min]');
    end
    
    subplot(3,2,5);
    plot(time_vec_minutes,larss_filter(:,iter_display),time_vec_minutes,larss_filter(:,iter_display),'*');
    title(['Larsson filter. F = ' num2str(F(iter_display))],'FontWeight','bold');
    xlabel('Time [Min]');
    ylabel('Amplitude');
    
    if ~Ignore_Gaussian_Calculation
        subplot(3,2,4);
        hold on;
        h1 = plot(time_vec_minutes,Sim_Ct_gauss_kernel(:,iter_display,iter_display),'r');
        h2 = plot(time_vec_minutes,Sim_Ct_gauss_kernel(:,iter_display,iter_display),'r*');
        h3 = plot(time_vec_minutes,Sim_Ct_gauss_kernel_noise(:,iter_display,iter_display),'b');
        h4 = plot(time_vec_minutes,Sim_Ct_gauss_kernel_noise(:,iter_display,iter_display),'bd');
        title('Gaussian filter result','FontWeight','bold');
        xlabel('Time [Min]');
        hold off;
        legend([h2 h4],'After Gauss filter','After Gauss filter+noise');
    end
    
    subplot(3,2,6);
    hold on;
    h1 = plot(time_vec_minutes,Sim_Ct_larss_kernel(:,iter_display,iter_display),'r');
    h2 = plot(time_vec_minutes,Sim_Ct_larss_kernel(:,iter_display,iter_display),'r*');
    h3 = plot(time_vec_minutes,Sim_Ct_larss_kernel_noise(:,iter_display,iter_display),'b');
    h4 = plot(time_vec_minutes,Sim_Ct_larss_kernel_noise(:,iter_display,iter_display),'bd');
    title('Larsson filter result','FontWeight','bold');
    ylabel('[mM]');
    xlabel('Time [Min]');
    
    % Display Larsson's parameters
    if Ignore_Gaussian_Calculation
        annotation('textbox',...
            [0.6 0.67 0.3 0.33],...
            'String',{'Larsson Filter Params:','',['F   = ' num2str(F(iter_display)) '  [mL/100g/min]'],...
            ['Vb = ' num2str(Vb_larss(iter_display)) '  [mL/100g]'],...
            ['Ve = ' num2str(Ve_larss(iter_display)) '  [mL/100g]'],...
            ['Vd = ' num2str(Vd(iter_display)) '  [mL/100g]'],...
            ['E   = ' num2str(E(iter_display))]},...
            'FontSize',8,...
            'FontName','Arial',...
            'LineStyle','-',...
            'EdgeColor',[0 0 0],...
            'LineWidth',2,...
            'BackgroundColor',[1 1 1],...
            'Color',[0.0 0.0 0]);
    else
        annotation('textbox',...
            [0.6 0.67 0.3 0.33],...
            'String',{'Larsson Filter Params:','',['F   = ' num2str(F(iter_display)) '  [mL/100g/min]'],...
            ['Vb = ' num2str(Vb_larss(iter_display)) '  [mL/100g]'],...
            ['Ve = ' num2str(Ve_larss(iter_display)) '  [mL/100g]'],...
            ['Vd = ' num2str(Vd(iter_display)) '  [mL/100g]'],...
            ['E   = ' num2str(E(iter_display))],'',...
            'Gauss Filter Params:','',...
            ['Delay Time = ' num2str(t_d(iter_display)*60) ' [Sec]'],...
            ['Sigma      = ' num2str(sigma(iter_display)*60) ' [Sec]'],...
            ['Amplitude  = ' num2str(amplitude(iter_display))]},...
            'FontSize',8,...
            'FontName','Arial',...
            'LineStyle','-',...
            'EdgeColor',[0 0 0],...
            'LineWidth',2,...
            'BackgroundColor',[1 1 1],...
            'Color',[0.0 0.0 0]);
    end
    
    
    hold off;
    legend([h2 h4],'After Larss filter','After Larss filter+noise');
    
    % Print result to PDF
    [Sim_Struct.idx_fig] = Print2Pdf(fig_num, Sim_Struct.idx_fig, 'Sim_AIF_Filters_Output.png', './Run_Output/',...
        'Sim AIF, Filters and Cts', 'SimulationGraphs');
    
end

time_finish = toc;
if ~strcmp(Verbosity,'None')
    display(sprintf('-I- Filtering AIF took %.2f seconds to finish...',time_finish));
end

if strcmp(Verbosity,'Full')
    display('-I- Finished filtering AIFs...');
end

end