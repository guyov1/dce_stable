function [ idx_fig ] = Plot_Simulation_Results( Sim_Struct,Verbosity, idx_fig)

% Take from struct variables used in local function
plot_error_results_flag           = Sim_Struct.plot_error_results_flag;
iterate_SNR                       = Sim_Struct.iterate_SNR;
iterate_sec_interval              = Sim_Struct.iterate_sec_interval;
Ignore_Gaussian_Calculation       = Sim_Struct.Ignore_Gaussian_Calculation;
iterate_gaussian_sigma            = Sim_Struct.iterate_gaussian_sigma;
iterate_gaussian_time_delay       = Sim_Struct.iterate_gaussian_time_delay;
iterate_gaussian_amplitude        = Sim_Struct.iterate_gaussian_amplitude;
iterate_F_larsson                 = Sim_Struct.iterate_F_larsson;
iterate_Vb_larsson                = Sim_Struct.iterate_Vb_larsson;
iterate_E_larsson                 = Sim_Struct.iterate_E_larsson;
Vb_single                         = Sim_Struct.Vb_single;
E_single                          = Sim_Struct.E_single;
F_single                          = Sim_Struct.F_single;
iterate_AIF_delay                 = Sim_Struct.iterate_AIF_delay;
results                           = Sim_Struct.results;

font_size                         = 20; % For plots

if strcmp(Verbosity,'Full')
    display('-I- Ploting simulation results...');
end

if (plot_error_results_flag)
    
    SNR_vec                     = results(1,:);
    time_interval_vec           = results(2,:);
    real_sigma_vec              = results(3,:);
    est_sigma_vec               = results(4,:);
    error_percent_sigma         = results(5,:);
    std_sigma                   = results(6,:);
    
    % Gaussian Delay Time
    real_t_d_vec                = results(7,:);  real_t_d_vec_sec = real_t_d_vec .* 60;
    est_t_d_vec                 = results(8,:);  est_t_d_vec_sec  = real_t_d_vec .* 60;
    error_percent_t_d           = results(9,:);
    std_t_d                     = results(10,:); std_t_d_sec      = std_t_d      .* 60;
    
    % Amp.
    real_amp_vec                = results(11,:);
    estimated_amp_vec           = results(12,:);
    error_percent_amp           = results(13,:);
    std_amp                     = results(14,:);
    
    % F
    real_larsson_F_vec          = results(15,:);
    est_larsson_F_vec           = results(16,:);
    error_percent_F             = results(17,:);
    std_F                       = results(18,:);
    
    % AIF Delay time (Larsson)
    real_t_d_Larss_vec_sec      = results(19,:);
    est_t_d_Larss_vec_sec       = results(20,:);
    error_percent_t_d_Larss     = results(21,:);
    std_t_d_Larss_sec           = results(22,:);
    
    % Ki
    real_larsson_Ki_Patlak_vec  = results(23,:);
    est_larsson_Ki_Patlak_vec   = results(24,:);
    error_percent_Ki_Patlak     = results(25,:);
    std_Ki_Patlak               = results(26,:);
    real_larsson_Ki_2CXM_vec    = results(27,:);
    est_larsson_Ki_2CXM_vec     = results(28,:);
    error_percent_Ki_2CXM       = results(29,:);
    std_Ki_2CXM                 = results(30,:);
    
    % PS
    real_larsson_PS_vec         = results(31,:);
    est_larsson_PS_vec          = results(32,:);
    error_percent_PS            = results(33,:);
    std_PS                      = results(34,:);
    
    % Vb
    real_larsson_Vb_Patlak_vec  = results(35,:);
    est_larsson_Vb_Patlak_vec   = results(36,:);
    error_percent_Vb_Patlak     = results(37,:);
    std_Vb_Patlak               = results(38,:);
    real_larsson_Vb_2CXM_vec    = results(39,:);
    est_larsson_Vb_2CXM_vec     = results(40,:);
    error_percent_Vb_2CXM       = results(41,:);
    std_Vb_2CXM                 = results(42,:);
    
    % Vd
    real_larsson_Vd_2CXM_vec    = results(43,:);
    est_larsson_Vd_2CXM_vec     = results(44,:);
    error_percent_Vd_2CXM       = results(45,:);
    std_Vd_2CXM                 = results(46,:);
    real_larsson_Vd_Normal_vec  = results(47,:);
    est_larsson_Vd_Normal_vec   = results(48,:);
    error_percent_Vd_Normal     = results(49,:);
    std_Vd_Normal               = results(50,:);
    
    % MTT
    real_larsson_MTT_2CXM_vec   = results(51,:);
    est_larsson_MTT_2CXM_vec    = results(52,:);
    error_percent_MTT_2CXM      = results(53,:);
    std_MTT_2CXM                = results(54,:);
    real_larsson_MTT_Normal_vec = results(55,:);
    est_larsson_MTT_Normal_vec  = results(56,:);
    error_percent_MTT_Normal    = results(57,:);
    std_MTT_Normal              = results(58,:);
    
    % E
    real_larsson_E_vec          = results(59,:);
    est_larsson_E_vec           = results(60,:);
    error_percent_E             = results(61,:);
    std_E                       = results(62,:);
    
    % AIF delay with Larss filter using Gaussian de-convolution
    real_t_d_Larss_vec_sec_using_Gauss      = results(63,:);
    est_t_d_Larss_using_Gauss_vec_sec       = results(64,:);
    error_percent_t_d_Larss_using_Gauss     = results(65,:);
    std_t_d_Larss_sec_using_Gauss           = results(66,:);
    
    % Sourbron - F
    real_larsson_F_Sourbron_vec         = results(67,:);
    est_larsson_F_Sourbron_vec          = results(68,:);
    error_percent_Sourbron_F            = results(69,:);
    std_Sourbron_F                      = results(70,:);
    
    % Sourbron - Ki
    real_larsson_Ki_Sourbron_2CXM_vec    = results(71,:);
    est_larsson_Ki_Sourbron_2CXM_vec     = results(72,:);
    error_percent_Ki_Sourbron_2CXM       = results(73,:);
    std_Ki_Sourbron_2CXM                 = results(74,:);
    
    % Sourbron - Vb
    real_larsson_Vb_Sourbron_2CXM_vec    = results(75,:);
    est_larsson_Vb_Sourbron_2CXM_vec     = results(76,:);
    error_percent_Vb_Sourbron_2CXM       = results(77,:);
    std_Vb_Sourbron_2CXM                 = results(78,:);
    
    % Sourbron - Ve
    real_larsson_Ve_Sourbron_2CXM_vec    = results(79,:);
    est_larsson_Ve_Sourbron_2CXM_vec     = results(80,:);
    error_percent_Ve_Sourbron_2CXM       = results(81,:);
    std_Ve_Sourbron_2CXM                 = results(82,:);
    
    
    if (iterate_SNR)
        
        if ~Ignore_Gaussian_Calculation
            
            fig_num = figure;
            
            subplot(3,1,1);
            plot(SNR_vec,error_percent_sigma);
            title('est. \sigma vs. SNR','FontWeight','bold');
            xlabel('SNR');
            ylabel('Error percent');
            
            subplot(3,1,2);
            plot(SNR_vec,error_percent_t_d);
            title('est. time delay vs. SNR','FontWeight','bold');
            xlabel('SNR');
            ylabel('Error percent');
            
            subplot(3,1,3);
            plot(SNR_vec,error_percent_amp);
            title('est. amplitude vs. SNR','FontWeight','bold');
            xlabel('SNR');
            ylabel('Error percent');
            
            % Print result to PDF
            [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Error_Analysis_Gauss_SNR.png', './Run_Output/',...
                'est. Error Analysis - SNR', 'EstErrorAnalysis');
            
        end
        
        % F versus noise
        
        fig_num = figure;
        
        subplot(3,2,1);
        plot(SNR_vec,error_percent_F);
        title(['est. F vs. SNR. F = ' num2str(real_larsson_F_vec(1))],'FontWeight','bold');
        xlabel('SNR');
        ylabel('Error percent');
        
        subplot(3,2,2);
        plot(SNR_vec,error_percent_E);
        title(['est. E vs. SNR. E = ' num2str(real_larsson_E_vec(1))],'FontWeight','bold');
        xlabel('SNR');
        ylabel('Error percent');
        
        subplot(3,2,3);
        plot(SNR_vec,error_percent_Vb_2CXM);
        title(['est. Vb vs. SNR. Vb = ' num2str(real_larsson_Vb_2CXM_vec(1))],'FontWeight','bold');
        xlabel('SNR');
        ylabel('Error percent');
        
        subplot(3,2,4);
        plot(SNR_vec,error_percent_Vd_2CXM);
        title(['est. Vd vs. SNR. Vd = ' num2str(real_larsson_Vd_2CXM_vec(1))],'FontWeight','bold');
        xlabel('SNR');
        ylabel('Error percent');
        
        subplot(3,2,5);
        plot(SNR_vec,error_percent_MTT_2CXM	);
        title(['est. MTT vs. SNR. MTT = ' num2str(real_larsson_MTT_2CXM_vec(1))],'FontWeight','bold');
        xlabel('SNR');
        ylabel('Error percent');
        
        subplot(3,2,6);
        plot(SNR_vec,error_percent_t_d_Larss);
        title(['est. DelayLarss vs. SNR. td = ' num2str(real_t_d_Larss_vec_sec(1))],'FontWeight','bold');
        xlabel('SNR');
        ylabel('Error percent');
        
        % Print result to PDF
        [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Error_Analysis_Larss_SNR.png', './Run_Output/',...
            'est. Error Analysis - SNR', 'EstErrorAnalysis');
        
        
        
    end
    
    if (iterate_sec_interval)
        
        if (~Ignore_Gaussian_Calculation)
            fig_num = figure;
            
            subplot(3,1,1);
            plot(time_interval_vec,error_percent_sigma);
            title('est. \sigma vs. time interval','FontWeight','bold');
            xlabel('Time Interval [sec]');
            ylabel('Error percent');
            
            subplot(3,1,2);
            plot(time_interval_vec,error_percent_t_d);
            title('est. amplitude vs. time interval','FontWeight','bold');
            xlabel('Time Interval [sec]');
            ylabel('Error percent');
            
            subplot(3,1,3);
            plot(time_interval_vec,error_percent_amp);
            title('est. time delay vs. time interval','FontWeight','bold');
            xlabel('Time Interval [sec]');
            ylabel('Error percent');
            
            % Print result to PDF
            [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Error_Analysis_time_interval.png', './Run_Output/',...
                'est. Error Analysis - Time Interval', 'EstErrorAnalysis');
        end
        
        fig_num = figure;
        
        subplot(2,1,1);
        %plot(real_sigma_vec,error_percent_sigma);
        errorbar(time_interval_vec,error_percent_F,std_F);
        title_string = 'est. error vs. time interval';
        title(title_string,'FontWeight','bold');
        xlabel('Time Interval [sec]');
        ylabel('Error percent');
        
        subplot(2,1,2);
        errorbar(time_interval_vec,est_larsson_F_vec,std_F);
        title_string = 'est. Flow vs. time interval';
        title(title_string,'FontWeight','bold');
        xlabel('Time Interval [sec]');
        ylabel('est. Flow');
        
        % Print result to PDF
        [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Error_Analysis_Larsson_F_vs_Time_Interval.png', './Run_Output/',...
            'est. F Error Analysis vs. Time Interval', 'EstErrorAnalysis');
        
        fig_num = figure;
        
        subplot(2,1,1);
        %plot(real_sigma_vec,error_percent_sigma);
        errorbar(time_interval_vec,error_percent_Ki_2CXM,std_Ki_2CXM);
        title_string = 'est. error vs. time interval';
        title(title_string,'FontWeight','bold');
        xlabel('Time Interval [sec]');
        ylabel('Error percent');
        
        subplot(2,1,2);
        errorbar(time_interval_vec,est_larsson_Ki_2CXM_vec,std_Ki_2CXM);
        title_string = 'est. Ki 2CXM vs. time interval';
        title(title_string,'FontWeight','bold');
        xlabel('Time Interval [sec]');
        ylabel('est. Ki 2CXM');
        
        % Print result to PDF
        [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Error_Analysis_Larsson_Ki_2CXM_vs_Time_Interval.png', './Run_Output/',...
            'est. Ki 2XCM Error Analysis vs. Time Interval', 'EstErrorAnalysis');
        
        
        fig_num = figure;
        
        subplot(2,1,1);
        %plot(real_sigma_vec,error_percent_sigma);
        errorbar(time_interval_vec,error_percent_MTT_2CXM,std_MTT_2CXM);
        title_string = 'est. error vs. time interval';
        title(title_string,'FontWeight','bold');
        xlabel('Time Interval [sec]');
        ylabel('Error percent');
        
        subplot(2,1,2);
        errorbar(time_interval_vec,est_larsson_MTT_2CXM_vec,std_MTT_2CXM);
        title_string = 'est. MTT 2CXM vs. time interval';
        title(title_string,'FontWeight','bold');
        xlabel('Time Interval [sec]');
        ylabel('est. MTT 2CXM');
        
        % Print result to PDF
        [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Error_Analysis_Larsson_MTT_2CXM_vs_Time_Interval.png', './Run_Output/',...
            'est. MTT 2XCM Error Analysis vs. Time Interval', 'EstErrorAnalysis');
        
    end
    
    % Display gaussian results in case we did not ignore it
    if (~Ignore_Gaussian_Calculation)
        
        if (iterate_gaussian_sigma)
            
            fig_num = figure;
            
            %plot(real_sigma_vec,error_percent_sigma);
            errorbar(real_sigma_vec,error_percent_sigma,std_sigma);
            title('est. error vs. original \sigma','FontWeight','bold');
            xlabel('Original sigma [sec]');
            ylabel('Error percent');
            
            % Print result to PDF
            [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Error_Analysis_Gaussian_sigma.png', './Run_Output/',...
                'est. Error Analysis - Original Gaussian Sigma Values Change', 'EstErrorAnalysis');
            
        end
        
        if (iterate_gaussian_time_delay)
            
            fig_num = figure;
            %plot(real_t_d_vec_sec,error_percent_t_d);
            errorbar(real_t_d_vec_sec,error_percent_t_d,std_t_d_sec);
            title('est. error vs. original time delay','FontWeight','bold');
            xlabel('Original time delay value [sec]');
            ylabel('Error percent');
            
            % Printing image to PDF
            gprint(fig_num,'Run_Output/Est_Error_Analysis_Gaussian_time_delay.png');
            
            idx_fig = idx_fig + 1;
            idx_string = ['idx_' num2str(idx_fig,'%03i')];
            AddToLog('Run_Output/',idx_string,'\\subsection*{\\underline{est. Error Analysis - Original Gaussian Time Delay Values Change}}');
            idx_fig = idx_fig + 1;
            idx_string = ['idx_' num2str(idx_fig,'%03i')];
            AddToLog('Run_Output/',idx_string,'EstErrorAnalysis','Est_Error_Analysis_Gaussian_time_delay.png');
            
        end
        
        if (iterate_gaussian_amplitude)
            
            fig_num = figure;
            %plot(real_sigma_vec,error_percent_sigma);
            errorbar(real_amp_vec,error_percent_amp,std_amp);
            title('est. error vs. original Amplitude','FontWeight','bold');
            xlabel('Original amplitude');
            ylabel('Error percent');
            
            % Print result to PDF
            [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Error_Analysis_Gaussian_Amp.png', './Run_Output/',...
                'est. Error Analysis - Original Gaussian Amplitude Values Change', 'EstErrorAnalysis');
            
        end
        
    end
    
    if (iterate_F_larsson)
        
        fig_num = figure;

        subplot(2,1,1);
        
        %plot(real_sigma_vec,error_percent_sigma);
        
        h1 = errorbar(real_larsson_F_vec,error_percent_F,std_F,'k');
        
        hold on;
        %h2 = errorbar(real_larsson_F_vec,error_percent_Sourbron_F,std_Sourbron_F,'r');
        
        %last_est_larsson_F_vec           = results(16,:);
        %last_error_percent_F             = results(17,:);
        %last_std_F                       = results(18,:);
        %save('\\fmri-guy2\Dropbox\University\Msc\Thesis\General\Matlab Simulations\Flow Extraction\Latest Code\Tmp.mat','last_est_larsson_F_vec','last_error_percent_F','last_std_F');
        
        %load('\\fmri-guy2\Dropbox\University\Msc\Thesis\General\Matlab Simulations\Flow Extraction\Latest Code\Tmp.mat');
        %h2 = errorbar(real_larsson_F_vec,last_error_percent_F,last_std_F,'g');
        hold off;
        
        %set(h_legend,'FontSize',font_size);
        
        title_string = sprintf('est. error vs. original Flow. Vb=%d, E=%.2f',Vb_single,E_single,'FontSize',font_size,'FontWeight','bold');
        title(title_string,'FontWeight','bold');
        xlabel('Original Flow','FontSize',font_size,'FontWeight','bold');
        ylabel('Error percent','FontSize',font_size,'FontWeight','bold');
        set(gca,'fontsize',font_size);
        
        subplot(2,1,2);
        
        errorbar(real_larsson_F_vec,est_larsson_F_vec,std_F,'b');
        hold on;
        errorbar(real_larsson_F_vec,est_larsson_F_Sourbron_vec,std_Sourbron_F,'r');
        %h2 = errorbar(real_larsson_F_vec,last_est_larsson_F_vec,last_std_F,'g');
        hold off;
        
        hr = refline(1); % y=x reference line
        set(hr,'Color','k','LineStyle','--','LineWidth',1.5);
        axis equal;   
        xlim([0 150]);  %Sets x axis limits
        ylim([0 150]);   %Sets y axis limits
        

        title_string = sprintf('est. Flow vs. original Flow. Vb=%d, E=%.2f',Vb_single,E_single,'FontSize',font_size,'FontWeight','bold');
        title(title_string,'FontWeight','bold');
        xlabel('Original Flow','FontSize',font_size,'FontWeight','bold');
        ylabel('est. Flow','FontSize',font_size,'FontWeight','bold');
        
        set(gca,'fontsize',font_size);
        
        % Print result to PDF
        [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Error_Analysis_Larsson_F.png', './Run_Output/',...
            'est. Error Analysis - Original Larsson Flow Values Change', 'EstErrorAnalysis');
        
    end
    
    if (iterate_Vb_larsson)
        
        fig_num = figure;
        
        subplot(2,1,1);
        %plot(real_sigma_vec,error_percent_sigma);
        errorbar(real_larsson_Vb_Patlak_vec,error_percent_Vb_Patlak,std_Vb_Patlak);
        title_string = sprintf('est. error vs. original Vb. F=%d, E=%.2f',F_single,E_single);
        title(title_string,'FontWeight','bold');
        xlabel('Original Vb');
        ylabel('Error percent');
        
        subplot(2,1,2);
        errorbar(real_larsson_Vb_Patlak_vec,est_larsson_Vb_Patlak_vec,std_Vb_Patlak);
        hr = refline(1); % y=x reference line
        set(hr,'Color','k','LineStyle','--','LineWidth',1.5);
        
        title_string = sprintf('est. Patlak Vb vs. original Vb. F=%d, E=%.2f',F_single,E_single);
        title(title_string,'FontWeight','bold');
        xlabel('Original Vb');
        ylabel('est. Vb');
        
        % Print result to PDF
        [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Error_Analysis_Patlak_Vb.png', './Run_Output/',...
            'est. Error Analysis - Original Larsson Vb Values Change', 'EstErrorAnalysis');
        
        fig_num = figure;
        
        subplot(2,1,1);
        %plot(real_sigma_vec,error_percent_sigma);
        errorbar(real_larsson_Vb_2CXM_vec,error_percent_Vb_2CXM,std_Vb_2CXM);
        title_string = sprintf('est. error vs. original Vb. F=%d, E=%.2f',F_single,E_single);
        title(title_string,'FontWeight','bold');
        xlabel('Original Vb');
        ylabel('Error percent');
        
        subplot(2,1,2);
        errorbar(real_larsson_Vb_2CXM_vec,est_larsson_Vb_2CXM_vec,std_Vb_2CXM);
        hr = refline(1); % y=x reference line
        set(hr,'Color','k','LineStyle','--','LineWidth',1.5);
        
        title_string = sprintf('est. 2CXM Vb vs. original Vb. F=%d, E=%.2f',F_single,E_single);
        title(title_string,'FontWeight','bold');
        xlabel('Original Vb');
        ylabel('est. Vb');
        
        % Print result to PDF
        [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Error_Analysis_Larsson_Vb.png', './Run_Output/',...
            'est. Error Analysis - Original Larsson Vb Values Change', 'EstErrorAnalysis');
        
    end
    
    if (iterate_E_larsson)
        
        fig_num = figure;
        
        subplot(2,1,1);
        %plot(real_sigma_vec,error_percent_sigma);
        errorbar(real_larsson_E_vec,error_percent_E,std_E);
        title_string = sprintf('est. error vs. original E. F=%d, Vb=%.2f',F_single,Vb_single);
        title(title_string,'FontWeight','bold');
        xlabel('Original E');
        ylabel('Error percent');
        
        subplot(2,1,2);
        errorbar(real_larsson_E_vec,est_larsson_E_vec,std_E);
        hr = refline(1); % y=x reference line
        set(hr,'Color','k','LineStyle','--','LineWidth',1.5);
        
        title_string = sprintf('est. E vs. original E. F=%d, Vb=%.2f',F_single,Vb_single);
        title(title_string,'FontWeight','bold');
        xlabel('Original E');
        ylabel('est. E');
        
        % Print result to PDF
        [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Error_Analysis_Patlak_E.png', './Run_Output/',...
            'est. Error Analysis - Original Larsson E Values Change', 'EstErrorAnalysis');
        
    end
    
    if (iterate_AIF_delay)
        
        fig_num = figure;
        
        %         subplot(2,1,1);
        %         %plot(real_sigma_vec,error_percent_sigma);
        %         errorbar(real_t_d_Larss_vec_sec,error_percent_t_d_Larss,std_t_d_Larss_sec,'b');
        %         if ~Ignore_Gaussian_Calculation
        %             hold on;
        %             errorbar(real_t_d_Larss_vec_sec_using_Gauss,error_percent_t_d_Larss_using_Gauss,std_t_d_Larss_sec_using_Gauss,'r');
        %             hold off;
        %         end
        %         title_string = sprintf('est. error vs. original AIF Delay. F=%d, Vb=%.2f. (RED-Gaussian)',F_single,Vb_single);
        %         title(title_string,'FontWeight','bold');
        %         xlabel('Original AIF delay');
        %         ylabel('Error percent');
        
        %subplot(2,1,2);
        errorbar(real_t_d_Larss_vec_sec,est_t_d_Larss_vec_sec,std_t_d_Larss_sec,'b');
        if ~Ignore_Gaussian_Calculation
            hold on;
            errorbar(real_t_d_Larss_vec_sec_using_Gauss,est_t_d_Larss_using_Gauss_vec_sec,std_t_d_Larss_sec_using_Gauss,'r');
            hold off;
        end
        
        hr = refline(1); % y=x reference line
        set(hr,'Color','k','LineStyle','--','LineWidth',1.5);
        
        title_string = sprintf('est. Larss AIF delay vs. original AIF delay. F=%d, Vb=%.2f. (RED-Gaussian)',F_single,Vb_single);
        title(title_string,'FontWeight','bold');
        xlabel('Original AIF delay');
        ylabel('est. AIF delay');
        
        % Print result to PDF
        [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Error_Analysis_AIF_Delay.png', './Run_Output/',...
            'est. Error Analysis - Original AIF Delay Values Change', 'EstErrorAnalysis');
        
        
        fig_num = figure;
        
        %plot(real_sigma_vec,error_percent_sigma);
        errorbar(real_t_d_Larss_vec_sec,error_percent_F,std_F);
        hold on;
        
        errorbar(real_t_d_Larss_vec_sec,error_percent_Sourbron_F,std_Sourbron_F,'r');
        
        hold off;
        
        title_string = sprintf('est. F error vs. AIF Delay. Vb=%d, E=%.2f, F=%.2f',Vb_single,E_single,real_larsson_F_vec(1));
        title(title_string,'FontWeight','bold');
        xlabel('AIF delay [Sec]');
        ylabel('Error percent');
        
        % Print result to PDF
        [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Est_Error_Analysis_Larsson_F_Vs_AIF_Delay.png', './Run_Output/',...
            'est. Error Analysis - Est. Larsson Flow error Vs. AIF Delay', 'EstErrorAnalysis');
        
        
    end
    
    
    
end


% Update return values
Sim_Struct.idx_fig = idx_fig;

end