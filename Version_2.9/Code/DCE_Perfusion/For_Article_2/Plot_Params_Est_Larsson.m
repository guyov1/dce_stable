close all;

font_size           = 26;
font_legend         = 24;
Line_Width          = 3;

min_Ktrans_display  = 0;
max_Ktrans_display  = 80;
min_E_display       = 0;
max_E_display       = 0.5;
min_BAT_display     = 0;
max_BAT_display     = 20;

Legend_Just_for_F   = true; % Plot legend just for flow parameter

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

idx_fig             = 1;
create_PDF          = false;
error_type          = 'rel' ; % Choose 'abs' or 'rel'
std_or_sem          = 'sem'; % Choose standard deviation or standard error of the mean
param_list_cells    = {'Flow' 'Ktrans' 'Vb' 'Ve' 'E' 'BAT'};
param_units_cells   = {'[mL/100g/min]' '[mL/100g/min]' '[mL/100g/min]' '[mL/100g/min]' '[a.u]' '[sec]'};

% ------------------------------------- Different Analysis Techniques -------------------------------------------

% DataPath_1          = [Base_Path 'Results_1000_Iterations_No_Delay_BiExp_Correction_2_sec_interval_6_min_total_Spline_2nd.mat'];
% DataPath_2          = [Base_Path 'Results_1000_Iterations_No_Delay_No_Correction_2_sec_interval_6_min_total_Spline_2nd.mat'];
% DataPath_3          = [Base_Path 'Results_1000_Iterations_No_Delay_No_Correction_2_sec_interval_6_min_total_PCA_2nd.mat'];
% DataPath_4          = [Base_Path 'Results_1000_Iterations_No_Delay_No_Correction_2_sec_interval_6_min_total_Wiener.mat'];
% Legend_1            = 'ACoPeD';
% Legend_2            = 'Spline 2nd';
% Legend_3            = 'PCA 2nd';
% Legend_4            = 'Wiener';
% E_y_scale           = 0.9/max_E_display;
% BAT_y_scale         = 1;
% num_data            = 4;
% Delay_Methods_flag  = false;

% ------------------------------------- with Delay -------------------------------------------

%DataPath_1          = [Base_Path 'Results_1000_Iterations_With_Delay_Up_To_20_With_Correction_2_sec_interval_6_min_total.mat'];
%DataPath_1          = [Base_Path 'Results_1000_Iterations_With_Delay_0_to_20_resolution_0.1_BiExp_Correction_2_sec_interval_6_min_total_Spline_2nd.mat'];
DataPath_1          = [Base_Path 'Results_1000_Iterations_With_Delay_0_to_20_resolution_0.05_BiExp_Correction_2_sec_interval_6_min_total_Spline_2nd.mat'];
%DataPath_2          = [Base_Path 'Results_1000_Iterations_With_Delay_Up_To_20_Cyclic_Correction_2_sec_interval_6_min_total.mat'];
DataPath_2          = [Base_Path 'Results_1000_Iterations_With_Delay_0_to_20_resolution_0.1_Cyclic_Correction_2_sec_interval_6_min_total_Spline_2nd.mat'];
%DataPath_3          = [Base_Path 'Results_1000_Iterations_With_Delay_0_to_20_resolution_0.1_Simple_Correction_2_sec_interval_6_min_total_Spline_2nd.mat'];
DataPath_3          = [Base_Path 'Results_1000_Iterations_With_Delay_0_to_20_resolution_0.1_LQ_MODEL_Correction_2_sec_interval_6_min_total_Spline_2nd.mat'];
DataPath_4          = [Base_Path 'Results_1000_Iterations_With_Delay_Up_To_20_No_Correction_2_sec_interval_6_min_total.mat'];
Legend_1            = 'ACoPeD';
Legend_2            = 'Cyclic Deconv. MVT';
%Legend_3            = 'MVT';
%Legend_2            = 'MVT';
Legend_3            = 'LQ';
Legend_4            = 'No Correct';
E_y_scale           = 1/max_E_display;
BAT_y_scale         = 1.4;
num_data            = 4;
Delay_Methods_flag  = true;

% ------------------------------------- Temproal Res. -------------------------------------------

% DataPath_1          = 'Results_1000_Iterations_No_Delay_BiExp_Correction_2_sec_interval_6_min_total_Spline_2nd.mat';
% DataPath_2          = 'Results_1000_Iterations_No_Delay_BiExp_Correction_4_sec_interval_6_min_total_Spline_2nd.mat';
% DataPath_3          = 'Results_1000_Iterations_No_Delay_BiExp_Correction_6_sec_interval_6_min_total_Spline_2nd.mat';
% Legend_1            = '2 sec';
% Legend_2            = '4 sec';
% Legend_3            = '6 sec';
% E_y_scale           = 0.65/max_E_display;
% BAT_y_scale         = 1;
% num_data            = 3;
% Delay_Methods_flag  = false;

% %------------------------------------- Scan duration. -------------------------------------------
% DataPath_1          = 'Results_1000_Iterations_No_Delay_BiExp_Correction_2_sec_interval_20_min_total_Spline_2nd.mat';
% DataPath_2          = 'Results_1000_Iterations_No_Delay_BiExp_Correction_2_sec_interval_6_min_total_Spline_2nd.mat';
% Legend_1            = '20 min';
% Legend_2            = '6 min';
% E_y_scale           = 0.5/max_E_display;
% BAT_y_scale         = 1;
% num_data            = 2;
% Delay_Methods_flag  = false;



legendInfo          = cell(1,num_data);
legendInfo_error    = cell(1,num_data);

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
    %eval(['AIF_delay_low_' num2str(idx) ' = AIF_delay_low;     ']);
    eval(['F_max_' num2str(idx) ' = F_max;      ']);
    eval(['E_max_' num2str(idx) ' = E_max;      ']);
    eval(['Vb_max_' num2str(idx) ' = Vb_max;     ']);
    eval(['Ve_max_' num2str(idx) ' = Ve_max;     ']);
    %eval(['AIF_delay_max_' num2str(idx) ' = AIF_delay_max;     ']);
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
minRealE                         = 0.07;
maxRealE                         = 0.5;
numBins                          = 15;
minVec                           = [minRealF minRealKtrans minRealVb minRealVe minRealE];
maxVec                           = [maxRealF maxRealKtrans maxRealVb maxRealVe maxRealE];
%generationMinVal                 = [F_low_1 E_low_1*F_low_1 Vb_low_1 Ve_low_1 E_low_1 AIF_delay_low_1];
%generationMaxVal                 = [F_max_1 E_max_1*F_max_1 Vb_max_1 Ve_max_1 E_max_1 AIF_delay_max_1];
generationMinVal                 = [F_low_1 E_low_1*F_low_1 Vb_low_1 Ve_low_1 E_low_1 ];
generationMaxVal                 = [F_max_1 E_max_1*F_max_1 Vb_max_1 Ve_max_1 E_max_1 ];

iter_idx = 0;

% Go over each parameter
for param_list = param_list_cells
    
    iter_idx = iter_idx + 1;    
    fig_num                          = figure;
    
    param_unit =  param_units_cells{iter_idx};
    param_name = param_list{:};
    
    
    
    hold on;
    % Dealing with Delay
    if Delay_Methods_flag && ( strcmp(param_name,'BAT') || strcmp(param_name,'Delay'))
        
%         for idx = [1 2 4]
         for idx = 1:num_data
            
            eval(['[real_bins_' num2str(idx) ', mean_bins_' num2str(idx) ', std_bins_' num2str(idx) ', sem_bins_' num2str(idx) ', mean_abs_error_' num2str(idx) ...
                ', std_abs_error_' num2str(idx) ', sem_abs_error_' num2str(idx) ', mean_rel_error_' num2str(idx) ', std_rel_error_' num2str(idx) ...
                ', sem_rel_error_' num2str(idx)' '] = binsDivision(realVec_' num2str(idx)...
                ', estVec_' num2str(idx) ', numBins, minVec, maxVec, generationMinVal, generationMaxVal, param_name);']);
            
            
            eval(['h' num2str(idx) ' = errorbar(real_bins_' num2str(idx) ',mean_bins_' num2str(idx) ',' std_or_sem '_bins_' num2str(idx) ...
                ',graph_format_' num2str(idx) ',''LineWidth'',Line_Width);']);
            
            
            %h1                               = scatter(real_bins, mean_bins,'g*');
            
            eval(['legendInfo{' num2str(idx) '} = Legend_' num2str(idx) ';']);
            
            eval(['mean_err   = mean_' error_type '_error_' num2str(idx) ';']);
            eval(['std_or_sem_str = ' std_or_sem '_' error_type '_error_' num2str(idx) ';']);
            eval(['legendInfo_error{' num2str(idx) '} = [Legend_' num2str(idx) ' '' '  num2str(mean_err,'%.2f') '+- ' num2str(std_or_sem_str,'%.2f') '''];']);
            
            legend boxoff;
        end
    
        % Dealing with all parameters except Delay
    else
        
        for idx = 1:num_data
            
            eval(['[real_bins_' num2str(idx) ', mean_bins_' num2str(idx) ', std_bins_' num2str(idx) ', sem_bins_' num2str(idx) ', mean_abs_error_' num2str(idx) ...
                ', std_abs_error_' num2str(idx) ', sem_abs_error_' num2str(idx) ', mean_rel_error_' num2str(idx) ', std_rel_error_' num2str(idx) ...
                ', sem_rel_error_' num2str(idx)' '] = binsDivision(realVec_' num2str(idx)...
                ', estVec_' num2str(idx) ', numBins, minVec, maxVec, generationMinVal, generationMaxVal, param_name);']);
            
            
            eval(['h' num2str(idx) ' = errorbar(real_bins_' num2str(idx) ',mean_bins_' num2str(idx) ',' std_or_sem '_bins_' num2str(idx) ',graph_format_' num2str(idx) ',''LineWidth'',Line_Width);']);
            
            
            %h1                               = scatter(real_bins, mean_bins,'g*');
            
            
            eval(['legendInfo{' num2str(idx) '} = Legend_' num2str(idx) ';']);
            
            eval(['mean_err   = mean_' error_type '_error_' num2str(idx) ';']);
            eval(['std_or_sem_str = ' std_or_sem '_' error_type '_error_' num2str(idx) ';']);
            eval(['legendInfo_error{' num2str(idx) '} = [Legend_' num2str(idx) ' '' '  num2str(mean_err,'%.2f') '+- ' num2str(std_or_sem_str,'%.2f') '''];']);
            legend boxoff; 
        end
        
    end
    
    
    hold off;
    
    % Plot legend just for F if required
    if (~Legend_Just_for_F || strcmp(param_name,'Flow'))
        [hleg, hobj] = legend(legendInfo,'Location','NorthWest');
        legend boxoff;
        textobj = findobj(hobj, 'type', 'text');
        set(textobj, 'Interpreter', 'latex', 'fontsize', font_legend);
    end
    
    %hr = refline(1); % y=x reference line
    %set(hr,'Color','k','LineStyle','--','LineWidth',Line_Width);
    
    title_string = sprintf(['Est. Vs. True ' param_name]);
    title(title_string,'FontSize',font_size,'FontWeight','bold');
    
    xlabel(['True ' param_name ' ' param_unit],'FontSize',font_size,'FontWeight','bold');
    ylabel(['Estimated ' param_name ' ' param_unit],'FontSize',font_size,'FontWeight','bold');
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
        axis([min_E_display, max_E_display,min_E_display, E_y_scale*max_E_display]);
    elseif ( strcmp(param_name,'BAT') || strcmp(param_name,'Delay'))
        axis([min_BAT_display, max_BAT_display,min_BAT_display, BAT_y_scale*max_BAT_display]);
    end
    
    hr = refline(1); % y=x reference line
    set(hr,'Color','k','LineStyle','--','LineWidth',Line_Width);
    saveas(fig_num,['./Run_Output/' param_name '_Comparison.jpg'])
    
    % Error Statistics
    fig_num       = figure;
    % mean_abs_error, std_abs_error, mean_rel_error, std_rel_error
    % error_type    = 'abs' ; % Choose 'abs' or 'rel'
    min_error_val = Inf;
    min_error_idx = NaN;
    hold on;
    for idx = 1:num_data
        eval(['h' num2str(idx) ' = errorbar(idx,mean_' error_type '_error_' num2str(idx) ',' std_or_sem '_' error_type '_error_' num2str(idx) ...
        ',graph_format_' num2str(idx) ',''LineWidth'',Line_Width);']);
        
        eval(['if (mean_' error_type '_error_' num2str(idx) ' < min_error_val)  min_error_val = mean_' error_type '_error_' ...
               num2str(idx) '  ; min_error_idx = ' num2str(idx) ' ; end']);
    
    end
    hold off;
    [hleg, hobj] = legend(legendInfo_error,'Location','NorthWest');
    legend boxoff;
    if ~isnan( min_error_idx )
        eval(['best_error = ' 'mean_' error_type '_error_' num2str(min_error_idx) ';']);
        eval(['best_' std_or_sem '   = ' std_or_sem '_'  error_type '_error_' num2str(min_error_idx) ';']);
    else
        best_error = Inf;
        eval(['best_' std_or_sem '   = NaN;']); 
    end
    
    if strcmp(std_or_sem,'sem')
        best_std_or_sem = best_sem;
    else
        best_std_or_sem = best_std;
    end

    title_string = sprintf(['Error Stats - ' param_name '. Best: ' num2str(best_error, '%.2f') '+-' num2str(best_std_or_sem, '%.2f') ]);
    
    title(title_string,'FontSize',font_size,'FontWeight','bold');
    xlabel(['True ' param_name ' ' param_unit],'FontSize',font_size,'FontWeight','bold');
    ylabel(['Error Val ' param_name ' ' param_unit],'FontSize',font_size,'FontWeight','bold');
    
    saveas(fig_num,['./Run_Output/' param_name '_Error_Comparison.jpg'])
    
    
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


close all;