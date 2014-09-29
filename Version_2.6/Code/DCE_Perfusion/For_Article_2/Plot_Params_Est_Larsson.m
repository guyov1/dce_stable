close all;

font_size           = 20;
font_legend         = 12;
Line_Width          = 3;

min_Ktrans_display  = 0;
max_Ktrans_display  = 80;
min_E_display       = 0;
max_E_display       = 0.5;
min_BAT_display     = 0;
max_BAT_display     = 3;
BAT_y_scale         = 2;


% Options
% '-'  - Solid line (default)
% '--' - Dashed line
% ':'  - Dotted line
% '-.' - Dash-dot line
% '+'	Plus sign
% 'o'	Circle
% '*'	Asterisk
% '.'	Point
% 'x'	Cross
% 'square' or 's' Square
% 'diamond' or 'd'Diamond
% '^'	Upward-pointing triangle
% 'v'	Downward-pointing triangle
% '>'	Right-pointing triangle
% '<'	Left-pointing triangle
% 'pentagram' or 'p' Five-pointed star (pentagram)
% 'hexagram' or 'h'  Six-pointed star (hexagram)
% r	Red
% g	Green
% b	Blue
% c	Cyan
% m	Magenta
% y	Yellow
% k	Black
% w	White

graph_format_1      = ':gs';
graph_format_2      = '--b+';
graph_format_3      = '-.c*';
graph_format_4      = '-ro';
Base_Path           = './Old_Runs/';
DataPath_1          = [Base_Path 'Results_1000_Iterations_With_Delay_No_Correction_2_sec_interval_6_min_total_Spline_2nd.mat'];
DataPath_2          = [Base_Path 'Results_1000_Iterations_With_Delay_Cyclic_Correction_2_sec_interval_6_min_total_Spline_2nd.mat'];
DataPath_3          = [Base_Path 'Results_1000_Iterations_With_Delay_And_Correction_2_sec_interval_6_min_total_Spline_2nd.mat'];
DataPath_4          = [Base_Path 'Results_1000_Iterations_With_Delay_And_Correction_2_sec_interval_6_min_total_Spline_2nd.mat'];
%DataPath_1          = 'Results_1000_Iterations_No_Delay_No_Correction_2_sec_interval_6_min_total.mat';
%DataPath_2          = 'Results_1000_Iterations_No_Delay_No_Correction_4_sec_interval_6_min_total.mat';
%DataPath_3          = 'Results_1000_Iterations_No_Delay_No_Correction_6_sec_interval_6_min_total.mat';
Legend_1            = 'MVT';
Legend_2            = 'Cyclic Deconv. MVT';
Legend_3            = 'BiExp';
Legend_4            = 'BiExp';
idx_fig             = 1;
create_PDF          = false;

param_list_cells    = {'Flow' 'Ktrans' 'Vb' 'Ve' 'E' 'Delay'};

%% Load all data

for idx = 1:num_data
    eval(['load(DataPath_' num2str(idx) ');']);
    eval(['plot_error_results_flag_' num2str(idx) ' = plot_error_results_flag;    ']);
    eval(['iterate_SNR_' num2str(idx) ' = iterate_SNR;']);
    eval(['iterate_sec_interval_' num2str(idx) ' = iterate_sec_interval;       ']);
    eval(['Ignore_Gaussian_Calculation_' num2str(idx) ' = Ignore_Gaussian_Calculation;']);
    eval(['iterate_gaussian_sigma_' num2str(idx) ' = iterate_gaussian_sigma;     ']);
    eval(['iterate_gaussian_time_delay_' num2str(idx) ' = iterate_gaussian_time_delay;']);
    eval(['iterate_gaussian_amplitude_' num2str(idx) ' = iterate_gaussian_amplitude; ']);
    eval(['iterate_F_larsson_' num2str(idx) ' = iterate_F_larsson;          ']);
    eval(['iterate_Vb_larsson_' num2str(idx) ' = iterate_Vb_larsson;         ']);
    eval(['iterate_E_larsson_' num2str(idx) ' = iterate_E_larsson;          ']);
    eval(['Vb_single_' num2str(idx) ' = Vb_single;  ']);
    eval(['E_single_' num2str(idx) ' = E_single;   ']);
    eval(['F_single_' num2str(idx) ' = F_single;   ']);
    eval(['iterate_AIF_delay_' num2str(idx) ' = iterate_AIF_delay;          ']);
    eval(['iterate_uniformly_' num2str(idx) ' = iterate_uniformly;          ']);
    eval(['Filter_Est_Chosen_' num2str(idx) ' = Filter_Est_Chosen;          ']);
    eval(['F_low_' num2str(idx) ' = F_low;      ']);
    eval(['E_low_' num2str(idx) ' = E_low;      ']);
    eval(['Vb_low_' num2str(idx) ' = Vb_low;     ']);
    eval(['Ve_low_' num2str(idx) ' = Ve_low;     ']);
    eval(['F_max_' num2str(idx) ' = F_max;      ']);
    eval(['E_max_' num2str(idx) ' = E_max;      ']);
    eval(['Vb_max_' num2str(idx) ' = Vb_max;     ']);
    eval(['Ve_max_' num2str(idx) ' = Ve_max;     ']);
    eval(['results_' num2str(idx) ' = results;    ']);
    eval(['results = results_' num2str(idx) ' ;      ']);
    eval(['real_larsson_F_vec_' num2str(idx) ' = results(15,:);  ']);
    eval(['est_larsson_F_vec_' num2str(idx) ' = results(16,:);  ']);
    eval(['error_percent_F_' num2str(idx) ' = results(17,:);  ']);
    eval(['std_F_' num2str(idx) ' = results(18,:);  ']);
    eval(['real_t_d_Larss_vec_sec_' num2str(idx) ' = results(19,:);  ']);
    eval(['est_t_d_Larss_vec_sec_' num2str(idx) ' = results(20,:);  ']);
    eval(['error_percent_t_d_Larss_' num2str(idx) ' = results(21,:);  ']);
    eval(['std_t_d_Larss_sec_' num2str(idx) ' = results(22,:);  ']);
    eval(['real_larsson_Ktrans_2CXM_vec_' num2str(idx) ' = results(27,:);  ']);
    eval(['est_larsson_Ktrans_2CXM_vec_' num2str(idx) ' = results(28,:);  ']);
    eval(['error_percent_Ktrans_2CXM_' num2str(idx) ' = results(29,:);  ']);
    eval(['std_Ktrans_2CXM_' num2str(idx) ' = results(30,:);  ']);
    eval(['real_larsson_Vb_Patlak_vec_' num2str(idx) ' = results(35,:);  ']);
    eval(['est_larsson_Vb_Patlak_vec_' num2str(idx) ' = results(36,:);  ']);
    eval(['error_percent_Vb_Patlak_' num2str(idx) ' = results(37,:);  ']);
    eval(['std_Vb_Patlak_' num2str(idx) ' = results(38,:);  ']);
    eval(['real_larsson_Vb_2CXM_vec_' num2str(idx) ' = results(39,:);  ']);
    eval(['est_larsson_Vb_2CXM_vec_' num2str(idx) ' = results(40,:);  ']);
    eval(['error_percent_Vb_2CXM_' num2str(idx) ' = results(41,:);  ']);
    eval(['std_Vb_2CXM_' num2str(idx) ' = results(42,:);  ']);
    eval(['real_larsson_E_vec_' num2str(idx) ' = results(59,:);  ']);
    eval(['est_larsson_E_vec_' num2str(idx) ' = results(60,:);  ']);
    eval(['error_percent_E_' num2str(idx) ' = results(61,:);  ']);
    eval(['std_E_' num2str(idx) ' = results(62,:);  ']);
    eval(['real_larsson_Ve_2CXM_vec_' num2str(idx) ' = results(83,:);  ']);
    eval(['est_larsson_Ve_2CXM_vec_' num2str(idx) ' = results(84,:);  ']);
    eval(['error_percent_Ve_2CXM_' num2str(idx) ' = results(85,:);  ']);
    eval(['std_Ve_2CXM_' num2str(idx) ' = results(86,:);  ']);
    eval(['real_F_' num2str(idx) '                            = real_larsson_F_vec_' num2str(idx) ';  ']);
    eval(['est_F_' num2str(idx) '                             = est_larsson_F_vec_' num2str(idx) ';  ']);
    eval(['real_Ktrans_' num2str(idx) '                       = real_larsson_Ktrans_2CXM_vec_' num2str(idx) ';  ']);
    eval(['est_Ktrans_' num2str(idx) '                        = est_larsson_Ktrans_2CXM_vec_' num2str(idx) ';  ']);
    eval(['real_Vb_' num2str(idx) '                           = real_larsson_Vb_2CXM_vec_' num2str(idx) ';  ']);
    eval(['est_Vb_' num2str(idx) '                            = est_larsson_Vb_2CXM_vec_' num2str(idx) ';  ']);
    eval(['real_Ve_' num2str(idx) '                           = real_larsson_Ve_2CXM_vec_' num2str(idx) ';  ']);
    eval(['est_Ve_' num2str(idx) '                            = est_larsson_Ve_2CXM_vec_' num2str(idx) ';  ']);
    eval(['real_E_' num2str(idx) '                            = real_larsson_E_vec_' num2str(idx) ';  ']);
    eval(['est_E_' num2str(idx) '                             = est_larsson_E_vec_' num2str(idx) ';  ']);
    eval(['real_BAT_' num2str(idx) '                          = real_t_d_Larss_vec_sec_' num2str(idx) ';  ']);
    eval(['est_BAT_' num2str(idx) '                           = est_t_d_Larss_vec_sec_' num2str(idx) ';  ']);
    eval(['realVec_' num2str(idx) '                           = [real_F_' num2str(idx) ' ;real_Ktrans_' num2str(idx) ' ;real_Vb_' num2str(idx) ' ;real_Ve_' num2str(idx) '; real_E_' num2str(idx) '; real_BAT_' num2str(idx) '];  ']);
    eval(['estVec_' num2str(idx) '                            = [est_F_' num2str(idx) ' ;est_Ktrans_' num2str(idx) ' ;est_Vb_' num2str(idx) ' ;est_Ve_' num2str(idx) '; est_E_' num2str(idx) '; est_BAT_' num2str(idx) '];  ']);
    
end

%% Plot
% --------------------------------------------
% Scatter plots for parameters
% --------------------------------------------

% Boundaries for F, Ktrans, Vb, Ve
minRealF                         = 0;
maxRealF                         = Inf;
minRealKtrans                    = 0;
maxRealKtrans                    = Inf;
minRealVb                        = 0;
maxRealVb                        = Inf;
minRealVe                        = 0;
maxRealVe                        = Inf;
minRealE                         = 0.05;
maxRealE                         = 0.5;
numBins                          = 15;
minVec                           = [minRealF minRealKtrans minRealVb minRealVe minRealE];
maxVec                           = [maxRealF maxRealKtrans maxRealVb maxRealVe maxRealE];
generationMinVal                 = [F_low_1 E_low_1*F_low_1 Vb_low_1 Ve_low_1 E_low_1 ];
generationMaxVal                 = [F_max_1 E_max_1*F_max_1 Vb_max_1 Ve_max_1 E_max_1 ];


for param_list = param_list_cells
    
    fig_num                          = figure;
    
    param_name = param_list{:};
    
    hold on;
    for idx = 1:num_data
        eval(['[real_bins_' num2str(idx) ', mean_bins_' num2str(idx) ', std_bins_' num2str(idx) '] = binsDivision(realVec_' num2str(idx)...
              ', estVec_' num2str(idx) ', numBins, minVec, maxVec, generationMinVal, generationMaxVal, param_name);']);
        eval(['h' num2str(idx) ' = errorbar(real_bins_' num2str(idx) ',mean_bins_' num2str(idx) ',std_bins_' num2str(idx) ...
              ',graph_format_' num2str(idx) ',''LineWidth'',Line_Width);']);
        
       %h1                               = scatter(real_bins, mean_bins,'g*');
          
        eval(['legendInfo{' num2str(idx) '} = Legend_' num2str(idx) ';']);
    end
    hold off;
    
    [hleg, hobj] = legend(legendInfo,'Location','NorthWest');
    
    
    textobj = findobj(hobj, 'type', 'text');
    set(textobj, 'Interpreter', 'latex', 'fontsize', font_legend);
    
    
    hr = refline(1); % y=x reference line
    set(hr,'Color','k','LineStyle','--','LineWidth',Line_Width);
    
    title_string = sprintf(['Est. Vs. Original ' param_name]);
    title(title_string,'FontSize',font_size,'FontWeight','bold');
    xlabel(['Original ' param_name],'FontSize',font_size,'FontWeight','bold');
    ylabel(['Estimated ' param_name],'FontSize',font_size,'FontWeight','bold');
    %legend([h1 h2], Legend_1, Legend_2);
    % Scale X and Y the same
    limx = get(gca, 'XLim');
    limy = get(gca, 'YLim');
    limBoth = [min([limx limy]) max([limx limy])];
    set(gca, 'XLim', limBoth,'FontSize',font_size,'FontWeight','bold');
    set(gca, 'YLim', limBoth,'FontSize',font_size,'FontWeight','bold');
    
    if strcmp(param_name,'Ktrans')
        axis([min_Ktrans_display, max_Ktrans_display,min_Ktrans_display, max_Ktrans_display]);
    elseif strcmp(param_name,'E')
        axis([min_E_display, max_E_display,min_E_display, max_E_display]);
    elseif ( strcmp(param_name,'BAT') || strcmp(param_name,'Delay'))
        axis([min_BAT_display, max_BAT_display,min_BAT_display, BAT_y_scale*max_BAT_display]);
    end
    
    saveas(fig_num,['./Run_Output/' param_name '_Comparison.jpg'])
    
    if create_PDF
        % Print result to PDF
        [idx_fig] = Print2Pdf(fig_num, idx_fig, ['Est_' param_name '_Uniform.png'], './Run_Output/',...
            [param_name ' Estimation - Uniform Generation'], ['Est' param_name]);
        
    end
    
    
end

if create_PDF
    % Create PDF Report
    Local_Path = [pwd filesep 'Run_Output'];
    Log_Path   = [pwd filesep 'Run_Output' filesep 'Log.mat'];
    %MakeReport_func(Local_Path, LogFN);
    MakeReport_func(Local_Path, Log_Path);
    
end


