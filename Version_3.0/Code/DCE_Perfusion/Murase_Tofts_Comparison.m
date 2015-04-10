function [ Sim_Struct ] = Murase_Tofts_Comparison( Sim_Struct, Verbosity )


if ~strcmp(Verbosity,'None')
    display('-I- Starting Murase/Tofts Comparison...');
end

% Take from struct variables used in local function
min_interval                  = Sim_Struct.min_interval;
time_vec_minutes              = Sim_Struct.time_vec_minutes;
Sim_AIF_delayed_no_noise      = Sim_Struct.Sim_AIF_delayed_no_noise;
Sim_AIF_with_noise            = Sim_Struct.Sim_AIF_with_noise;
SNR_ratio                     = Sim_Struct.SNR_ratio;
One_Iteration_Murase_Tofts    = Sim_Struct.One_Iteration_Murase_Tofts;
Iterate_Murase_Tofts_num_iter = Sim_Struct.Iterate_Murase_Tofts_num_iter;
algorithm_options             = Sim_Struct.algorithm_options;
adjusted_larsson              = Sim_Struct.Adjusted_Larsson_Model;

% Check if only one iteration is wanted
if (One_Iteration_Murase_Tofts)
    Iterate_Murase_Tofts_num_iter = 1;
end

% AIF index to choose for iteration
AIF_idx = 1;

% Initiate parameter matrices
Murase_params    = zeros(3,Iterate_Murase_Tofts_num_iter);
NonLinear_params = zeros(3,Iterate_Murase_Tofts_num_iter);
Perfusion_params = zeros(4,Iterate_Murase_Tofts_num_iter);

% Initiate parameters for simulation
F_low_val         = 10;
F_max_val         = 150;
Vb_larss_Murase   = 6;                                % [mL/100g]     , they used 3,6,12,18
Ve_larss_Murase   = 2; % Must be smaller than Vtis
Hct               = 0.38;

% Initiate time calculation for first index
tic;
display(sprintf('-I- Starting simulation for %d voxels...',Iterate_Murase_Tofts_num_iter));

%parfor idx = 1 : Iterate_Murase_Tofts_num_iter
for idx = 1 : Iterate_Murase_Tofts_num_iter
    
    % Report after each 100 voxels
    if ( mod(idx,1000) == 0 )
        display(sprintf('Finished Regularization_Methods_Simulation for %d voxels...',idx));
        %remaining_voxels = num_voxels - j;
        %fprintf('Number of remaining voxels: %d .\n',remaining_voxels);
    end
    
    % Values for a single iteration
    F_Murase               = 120;
    E_Murase               = 0.0;                              % Extraction. They used 0,0.1,0.5,1
    
    if (~One_Iteration_Murase_Tofts)
        F_Murase               = F_low_val + (F_max_val-F_low_val)*rand;
    end
    
    %E_Murase               = 0.10;
    % Extraction. They used 0,0.1,0.5,1
    if (~One_Iteration_Murase_Tofts)
        E_Murase               = rand;                              % Extraction. They used 0,0.1,0.5,1
    end
    
    Ktrans_Murase              = E_Murase * F_Murase;
    %Ve_larss_Murase        = 0.10 * ( 1 - (Vb_larss_Murase/100) )*100; % Must be smaller than Vtis
    
    if adjusted_larsson
        IRF_larss_Murase       = Adjusted_Larsson_Filter(time_vec_minutes, F_Murase, Vb_larss_Murase, E_Murase, Ve_larss_Murase);  % No units
    else
        IRF_larss_Murase       = Larsson_Filter(time_vec_minutes, F_Murase, Vb_larss_Murase, E_Murase, Ve_larss_Murase, Hct);  % No units
    end
    
    
    larss_filter_Murase    = F_Murase*IRF_larss_Murase; % [mL/100g/min]
    
    % Filter the delayed AIF with Larsson's kernel
    if Sim_Struct.ignore_time_delta
        Sim_Ct_larss_Murase = filter(larss_filter_Murase,1,Sim_AIF_delayed_no_noise(:,AIF_idx,AIF_idx));
    else
        Sim_Ct_larss_Murase = filter(larss_filter_Murase*min_interval(AIF_idx),1,Sim_AIF_delayed_no_noise(:,AIF_idx,AIF_idx));
    end
    
    % Adding noise to simulated Ct(t) from both kernels
    noise_sigma_larss_Murase   = mean(Sim_Ct_larss_Murase)/SNR_ratio(AIF_idx);
    noise_to_add_larss_Murase  = noise_sigma_larss_Murase * randn(size(Sim_Ct_larss_Murase));
    Sim_Ct_larss_Murase_noise = Sim_Ct_larss_Murase + noise_to_add_larss_Murase;
    % Zero negative values (non realistic)
    Sim_Ct_larss_Murase_noise(Sim_Ct_larss_Murase_noise<0) = 0;
    
    % Create vectors of A matrix
    A_1 =  cumtrapz(time_vec_minutes,Sim_AIF_with_noise(:,AIF_idx,AIF_idx));
    A_2 = -cumtrapz(time_vec_minutes,Sim_Ct_larss_Murase_noise);
    A_3 =  Sim_AIF_with_noise(:,AIF_idx,AIF_idx);
    A   =  [A_1' A_2' A_3'];
    
    % Create c vector
    C_vec = Sim_Ct_larss_Murase_noise';
    B     = A \ C_vec;
    
    K_trans_Tofts_Murase = B(1) - B(2)*B(3);
    K_ep_Tofts_Murase    = B(2);
    Vp_Tofts_Murase      = B(3);
    
    
    % Non-Linear fit
    % The analytic funcational of Tofts model
    Tofts_function    = @(x,t) Tofts_model( t, min_interval(AIF_idx), Sim_AIF_with_noise(:,AIF_idx,AIF_idx), x(1), x(2), x(3));
    % Set Murase as initial start point for non-linear search
    Init_Guess_Tofts  = [K_trans_Tofts_Murase   K_ep_Tofts_Murase   Vp_Tofts_Murase]; %Ktrans, Kep ,Vp
    LowerBound_Tofts  = [0   0   0 ];
    UpperBound_Tofts  = [100 100 100];
    
    [est_params_Tofts,residue_norm_Tofts,residual_Tofts,exitflag_Tofts,algo_info_Tofts] = ...
        lsqcurvefit(Tofts_function,Init_Guess_Tofts,time_vec_minutes',Sim_Ct_larss_Murase_noise,...
        LowerBound_Tofts,UpperBound_Tofts,algorithm_options);
    
    K_trans_Tofts_Murase_NonLinear = est_params_Tofts(1);
    K_ep_Tofts_Murase_NonLinear    = est_params_Tofts(2);
    Vp_Tofts_Murase_NonLinear      = est_params_Tofts(3);
    
    
    % Put params in vectors
    %Murase_params(1,idx)    = K_trans_Tofts_Murase;
    %Murase_params(2,idx)    = K_ep_Tofts_Murase;
    %Murase_params(3,idx)    = Vp_Tofts_Murase;
    Murase_params(:,idx)    = [K_trans_Tofts_Murase K_ep_Tofts_Murase Vp_Tofts_Murase];
    
    %NonLinear_params(1,idx) = K_trans_Tofts_Murase_NonLinear;
    %NonLinear_params(2,idx) = K_ep_Tofts_Murase_NonLinear;
    %NonLinear_params(3,idx) = Vp_Tofts_Murase_NonLinear;
    NonLinear_params(:,idx) = [K_trans_Tofts_Murase_NonLinear K_ep_Tofts_Murase_NonLinear Vp_Tofts_Murase_NonLinear];
    
    %Perfusion_params(1,idx) = F_Murase;
    %Perfusion_params(2,idx) = Vb_larss_Murase;
    %Perfusion_params(3,idx) = E_Murase;
    %Perfusion_params(4,idx) = Ve_larss_Murase;
    Perfusion_params(:,idx) = [F_Murase Vb_larss_Murase E_Murase Ve_larss_Murase];
    
    
end

display(sprintf('Finished simulation for %d voxels...',Iterate_Murase_Tofts_num_iter));
time_finish = toc;
display(sprintf('Took %.2f seconds to finish...',time_finish));
tic;

Real_F_vec     = Perfusion_params(1,:);
Real_E_vec     = Perfusion_params(3,:);

Est_Ktrans_vec = Murase_params(1,:);
Est_Kep_vec    = Murase_params(2,:);
Est_Vp_vec     = Murase_params(3,:);

Est_Ktrans_NonLinear_vec = NonLinear_params(1,:);
Est_Kep_NonLinear_vec    = NonLinear_params(2,:);
Est_Vp_NonLinear_vec     = NonLinear_params(3,:);

if (One_Iteration_Murase_Tofts) % Plot
    
    % Murase fit
    AIF_part                       = Est_Vp_vec*Sim_AIF_with_noise(:,AIF_idx,AIF_idx);
    if Sim_Struct.ignore_time_delta
        Kep_Filter_Part                = filter(Est_Ktrans_vec*exp(-Est_Kep_vec*time_vec_minutes),1,Sim_AIF_with_noise(:,AIF_idx,AIF_idx));
    else
        Kep_Filter_Part                = filter(Est_Ktrans_vec*exp(-Est_Kep_vec*time_vec_minutes)*min_interval(AIF_idx),1,Sim_AIF_with_noise(:,AIF_idx,AIF_idx));
    end
    
    Sum_Result                     = AIF_part + Kep_Filter_Part;
    Squared_Error_vs_time          = (Sim_Ct_larss_Murase_noise - Sum_Result).^2;
    Squared_Error_vs_time_Kep_Only = (Sim_Ct_larss_Murase_noise - Kep_Filter_Part).^2;
    
    % Non linear fit
    AIF_part_NonLinear                       = Est_Vp_NonLinear_vec*Sim_AIF_with_noise(:,AIF_idx,AIF_idx);
    
    if Sim_Struct.ignore_time_delta
        Kep_Filter_Part_NonLinear                = filter(Est_Ktrans_NonLinear_vec*exp(-Est_Kep_NonLinear_vec*time_vec_minutes),1,Sim_AIF_with_noise(:,AIF_idx,AIF_idx));
    else
        Kep_Filter_Part_NonLinear                = filter(Est_Ktrans_NonLinear_vec*exp(-Est_Kep_NonLinear_vec*time_vec_minutes)*min_interval(AIF_idx),1,Sim_AIF_with_noise(:,AIF_idx,AIF_idx));
    end
    
    Sum_Result_NonLinear                     = AIF_part_NonLinear + Kep_Filter_Part_NonLinear;
    Squared_Error_vs_time_NonLinear          = (Sim_Ct_larss_Murase_noise - Sum_Result_NonLinear).^2;
    Squared_Error_vs_time_Kep_Only_NonLinear = (Sim_Ct_larss_Murase_noise - Kep_Filter_Part_NonLinear).^2;
    Squared_Error_vs_time_AIF_Only_NonLinear = (Sim_Ct_larss_Murase_noise - AIF_part_NonLinear).^2;
    
    
    fig_num = figure;
    subplot(3,3,1);
    plot(time_vec_minutes,Sim_AIF_with_noise(:,AIF_idx,AIF_idx),'k');
    title('Arterial Input Function');
    ylabel('C_a(t) [mM]');
    xlabel('Time [Min]');
    % Murase
    subplot(3,3,4);
    hold on;
    h1 = plot(time_vec_minutes,larss_filter_Murase,'LineWidth',5,'Color','k');
    h2 = plot(time_vec_minutes,Est_Ktrans_vec*exp(-Est_Kep_vec*time_vec_minutes),'LineWidth',5,'LineStyle',':','Color','b');
    h3 = stem(0,Est_Vp_vec,'fill','LineWidth',5,'Color','r');
    hold off;
    title('Larsson vs Tofts IRF');
    xlabel('Time [Min]');
    ylabel('IRF');
    legend([h1 h2 h3],'Larsson IRF','Tofts IRF - permeability','Tofts IRF - AIF')
    subplot(3,3,5);
    hold on;
    h1 = plot(time_vec_minutes,Sim_Ct_larss_Murase_noise,'LineWidth',6,'Color','k');
    h2 = plot(time_vec_minutes,AIF_part,'LineWidth',1,'Marker','+','Color','r');
    h3 = plot(time_vec_minutes,Kep_Filter_Part,'LineWidth',1,'LineStyle','o','Color','g');
    h4 = plot(time_vec_minutes,Sum_Result,'LineWidth',2,'Color','b');
    %h5 = plot(time_vec_minutes,Sum_Result_NonLinear,'LineWidth',2,'Color','c');
    hold off;
    title('Concentration Time Curve');
    legend([h1 h2 h3 h4],'CTC','Tofts AIF','Tofts permeability','Tofts fit')
    xlabel('Time [Min]');
    ylabel('C_t(t) [mM]');
    subplot(3,3,6);
    hold on;
    h1 = plot(time_vec_minutes,Squared_Error_vs_time,'k');
    h2 = plot(time_vec_minutes, cumtrapz(time_vec_minutes,Squared_Error_vs_time),'k+');
    h3 = plot(time_vec_minutes,Squared_Error_vs_time_Kep_Only,'b--');
    h4 = plot(time_vec_minutes, cumtrapz(time_vec_minutes,Squared_Error_vs_time_Kep_Only),'bo');
    h_legend = legend([h1 h2 h3 h4],'Tofts Error(ti)','Cummulative Tofts Error',...
        'Tofts Permeability Error(ti)','Cummulative Tofts Permeability Error');
    set(h_legend,'FontSize',8);
    title('Squared Error');
    xlabel('Time [Min]');
    ylabel('Error');
    % Non Linear Fit
    subplot(3,3,7);
    hold on;
    h1 = plot(time_vec_minutes,larss_filter_Murase,'LineWidth',5,'Color','k');
    h2 = plot(time_vec_minutes,Est_Ktrans_NonLinear_vec*exp(-Est_Kep_NonLinear_vec*time_vec_minutes),'LineWidth',5,'LineStyle',':','Color','b');
    h3 = stem(0,Est_Vp_NonLinear_vec,'fill','LineWidth',5,'Color','r');
    hold off;
    title('Larsson vs Tofts IRF - Non Linear');
    xlabel('Time [Min]');
    ylabel('IRF');
    legend([h1 h2 h3],'Larsson IRF','Tofts IRF - permeability','Tofts IRF - AIF');
    subplot(3,3,8);
    hold on;
    h1 = plot(time_vec_minutes,Sim_Ct_larss_Murase_noise,'LineWidth',6,'Color','k');
    h2 = plot(time_vec_minutes,AIF_part_NonLinear,'LineWidth',1,'LineStyle','+','Color','r');
    h3 = plot(time_vec_minutes,Kep_Filter_Part_NonLinear,'LineWidth',1,'LineStyle','o','Color','g');
    h4 = plot(time_vec_minutes,Sum_Result_NonLinear,'LineWidth',2,'Color','b');
    %h5 = plot(time_vec_minutes,Sum_Result_NonLinear,'LineWidth',2,'Color','c');
    hold off;
    title('Concentration Time Curve');
    legend([h1 h2 h3 h4],'CTC','Tofts AIF','Tofts permeability','Tofts fit')
    xlabel('Time [Min]');
    ylabel('C_t(t) [mM]');
    subplot(3,3,9);
    hold on;
    h1 = plot(time_vec_minutes,Squared_Error_vs_time_NonLinear,'k');
    h2 = plot(time_vec_minutes, cumtrapz(time_vec_minutes,Squared_Error_vs_time_NonLinear),'k+');
    h3 = plot(time_vec_minutes,Squared_Error_vs_time_AIF_Only_NonLinear,'b--');
    h4 = plot(time_vec_minutes, cumtrapz(time_vec_minutes,Squared_Error_vs_time_AIF_Only_NonLinear),'bo');
    %h3 = plot(time_vec_minutes,Squared_Error_vs_time_Kep_Only_NonLinear,'b--');
    %h4 = plot(time_vec_minutes, cumtrapz(time_vec_minutes,Squared_Error_vs_time_Kep_Only_NonLinear),'bo');
    
    h_legend = legend([h1 h2 h3 h4],'Tofts Error(ti)','Cummulative Tofts Error',...
        'Tofts AIF Error(ti)','Cummulative Tofts AIF Error');
    set(h_legend,'FontSize',8);
    title('Squared Error');
    xlabel('Time [Min]');
    ylabel('Error');
    
    
    % Display Larsson's parameters
    annotation('textbox',...
        [0.7 0.75 0.2 0.23],...
        'String',{...
        ['Temporal res. = ' num2str(sec_interval) '  [sec]'],...
        ['Flow          = ' num2str(F_Murase) '  [mL/100g/min]'],...
        ['Extraction    = ' num2str(E_Murase,'%.2f')],...
        ['Vb            = ' num2str(Vb_larss_Murase) '  [mL/100g]'],...
        ['Ve            = ' num2str(Ve_larss_Murase,'%.2f') '  [mL/100g]']},...
        'FontSize',8,...
        'FontName','Arial',...
        'LineStyle','-',...
        'EdgeColor',[0 0 0],...
        'LineWidth',2,...
        'BackgroundColor',[1 1 1],...
        'Color',[0.0 0.0 0]);
    
    % Print result to PDF
    [Sim_Struct.idx_fig] = Print2Pdf(fig_num, Sim_Struct.idx_fig, 'Abstract_Analysis.png', './Run_Output/',...
        'Analysis when E=0', 'AnalysisNoPermeability');
    
    % More specified plot for abstract
    num_rows_fig = 1;
    num_cols_fig = 2;
    subplot_idx  = 1;
    
    % Display only 2 mintues (so graphs will be thinner)
    Display_minute = 2;
    to_keep_ration = Display_minute / total_sim_time_min;
    total_samples  = max(size(time_vec_minutes));
    disp_samples   = round(to_keep_ration*total_samples);
    
    fig_num      = figure;
    font_size    = 25;
    line_size_1  = 10;
    line_size_2  = 10;
    
    %                 subplot(num_rows_fig,num_cols_fig,subplot_idx);
    %                 subplot_idx = subplot_idx + 1;
    %                 plot(time_vec_minutes,Sim_AIF_with_noise(:,AIF_idx,AIF_idx),'k');
    %                 title('Arterial Input Function');
    %                 ylabel('C_a(t) [mM]');
    %                 xlabel('Time [Min]');
    %
    subplot(num_rows_fig,num_cols_fig,subplot_idx);
    subplot_idx = subplot_idx + 1;
    set(gca,'fontsize',10,'FontWeight','bold');
    hold on;
    h1 = plot(time_vec_minutes(1:disp_samples),larss_filter_Murase(1:disp_samples),'LineWidth',line_size_1,'Color','k');
    h2 = plot(time_vec_minutes(1:disp_samples),Est_Ktrans_NonLinear_vec*exp(-Est_Kep_NonLinear_vec*time_vec_minutes(1:disp_samples)),'LineWidth',line_size_1/2,'LineStyle','o','Color','b');
    h3 = stem(0,Est_Vp_NonLinear_vec,'fill','LineWidth',5,'Color','r');
    hold off;
    title('Larsson vs Tofts IRF','FontSize',font_size,'FontWeight','bold');
    xlabel('Time [Min]','FontSize',font_size,'FontWeight','bold');
    ylabel('IRF(t)','FontSize',font_size,'FontWeight','bold');
    set(gca,'fontsize',font_size);
    h_legend = legend([h1 h2 h3],'Larsson IRF','Tofts IRF - permeability','Tofts IRF - AIF');
    set(h_legend,'FontSize',font_size);
    
    subplot(num_rows_fig,num_cols_fig,subplot_idx);
    subplot_idx = subplot_idx + 1;
    set(gca,'fontsize',font_size,'FontWeight','bold');
    hold on;
    AIF_to_plot = Sim_AIF_with_noise(:,AIF_idx,AIF_idx);
    h1 = plot(time_vec_minutes(1:disp_samples),AIF_to_plot(1:disp_samples),'LineWidth',line_size_2,'LineStyle','--','Color','c');
    h2 = plot(time_vec_minutes(1:disp_samples),Sim_Ct_larss_Murase_noise(1:disp_samples),'LineWidth',line_size_2,'Color','k');
    h3 = plot(time_vec_minutes(1:disp_samples),AIF_part_NonLinear(1:disp_samples),'LineWidth',line_size_2/4,'LineStyle','+','Color','r');
    h4 = plot(time_vec_minutes(1:disp_samples),Kep_Filter_Part_NonLinear(1:disp_samples),'LineWidth',line_size_2/4,'LineStyle','o','Color','g');
    h5 = plot(time_vec_minutes(1:disp_samples),Sum_Result_NonLinear(1:disp_samples),'LineWidth',line_size_2/2,'LineStyle','-.','Color','m');
    %h5 = plot(time_vec_minutes,Sum_Result_NonLinear,'LineWidth',2,'Color','c');
    hold off;
    title('Concentration Time Curve','FontSize',font_size,'FontWeight','bold');
    h_legend = legend([h1 h2 h3 h4 h5],'AIF','CTC','Tofts AIF','Tofts permeability','Tofts fit');
    set(h_legend,'FontSize',font_size);
    xlabel('Time [Min]','FontSize',font_size,'FontWeight','bold');
    ylabel('C_t(t) [mM]','FontSize',font_size,'FontWeight','bold');
    set(gca,'fontsize',font_size);
    
    % Save figure as jpeg
    saveas(fig_num,'\\fmri-guy2\Dropbox\University\Msc\Thesis\ISMRM 2014\Mine\Larsson_Vs_Tofts.jpg');
    
    %                 subplot(num_rows_fig,num_cols_fig,subplot_idx);
    %                 subplot_idx = subplot_idx + 1;
    %                 hold on;
    %
    %                 h1 = plot(time_vec_minutes,Squared_Error_vs_time_NonLinear,'k');
    %                 h2 = plot(time_vec_minutes, cumtrapz(time_vec_minutes,Squared_Error_vs_time_NonLinear),'k+');
    %                 h3 = plot(time_vec_minutes,Squared_Error_vs_time_AIF_Only_NonLinear,'b--');
    %                 h4 = plot(time_vec_minutes, cumtrapz(time_vec_minutes,Squared_Error_vs_time_AIF_Only_NonLinear),'bo');
    %                 %h3 = plot(time_vec_minutes,Squared_Error_vs_time_Kep_Only_NonLinear,'b--');
    %                 %h4 = plot(time_vec_minutes, cumtrapz(time_vec_minutes,Squared_Error_vs_time_Kep_Only_NonLinear),'bo');
    %
    %                 h_legend = legend([h1 h2 h3 h4],'Tofts Error(ti)','Cummulative Tofts Error',...
    %                     'Tofts AIF Error(ti)','Cummulative Tofts AIF Error');
    %
    %                 set(h_legend,'FontSize',8);
    %                 title('Squared Error - Non Linear');
    %                 xlabel('Time [Min]');
    %                 ylabel('Error(t)');
    
    % Print result to PDF
    [Sim_Struct.idx_fig] = Print2Pdf(fig_num, Sim_Struct.idx_fig, 'Abstract_Analysis_Smaller.png', './Run_Output/',...
        'Analysis when E=0 for abstract', 'AnalysisNoPermeability4Abstract');
    
else
    
    
    
    % --------------   Ktrans (vs. E and F)    --------------
    fig_num = figure;
    
    tri = delaunay(Real_F_vec, Real_E_vec);
    trisurf(tri,Real_F_vec,Real_E_vec,Est_Ktrans_vec);
    % Make the plot look nicer
    shading interp;
    %xlabel('F','Color',[0 0 1]);
    hold all;
    sorted_F = sort(Real_F_vec);
    z_straight_line = 18 + ( (30-10)/(sorted_F(end) - sorted_F(1) ) )*sorted_F ;
    h = plot3(sorted_F,repmat(0.06,1,Iterate_Murase_Tofts_num_iter),z_straight_line,'--k','LineWidth',3);
    %plot3(sorted_F,0.12,32,'*k','LineWidth',5);
    hold off;
    
    xlabel('F[mL/100mL/min]','fontsize',25,'Color','k','FontWeight','bold');
    ylabel('E','fontsize',25,'Color','k','FontWeight','bold');
    zlabel('Ktrans[mL/100mL/min]','fontsize',25,'Color','k','FontWeight','bold');
    set(gca,'fontsize',23,'FontWeight','bold');
    title(sprintf('Time resolution = %d sec. Vb = %.2f, Ve = %.2f.',sec_interval,...
        Vb_larss_Murase,Ve_larss_Murase),'fontsize',25,'FontWeight','bold');
    legend(h,'E = 0.06','fontsize',25);
    h_color = colorbar;
    set(h_color,'fontsize',23,'FontWeight','bold');
    
    % Print result to PDF
    [Sim_Struct.idx_fig] = Print2Pdf(fig_num, Sim_Struct.idx_fig, 'Ktrans_vs_E_F.png', './Run_Output/',...
        'Ktrans vs. E and F', 'KtransVsEandF');
    
    %  --------------  Kep (vs. E and F)    --------------
    fig_num = figure;
    tri = delaunay(Real_F_vec, Real_E_vec);
    trisurf(tri,Real_F_vec,Real_E_vec,Est_Kep_vec);
    % Make the plot look nicer
    shading interp;
    %xlabel('F','Color',[0 0 1]);
    
    xlabel('F [mL/100mL/min]','Color',[0 0 1]);
    ylabel('E','Color',[0 0 1]);
    zlabel('Kep [mL/100mL/min]','Color',[0 0 1]);
    title(sprintf('Time resolution = %d sec. Vb = %.2f, Ve = %.2f.',sec_interval,Vb_larss_Murase,Ve_larss_Murase));
    
    % Print result to PDF
    [Sim_Struct.idx_fig] = Print2Pdf(fig_num, Sim_Struct.idx_fig, 'Kep_vs_E_F.png', './Run_Output/',...
        'Kep vs. E and F', 'KepVsEandF');
    
    %   -------------- Vp (vs. E and F)    --------------
    fig_num = figure;
    tri = delaunay(Real_F_vec, Real_E_vec);
    trisurf(tri,Real_F_vec,Real_E_vec,Est_Vp_vec);
    % Make the plot look nicer
    shading interp;
    %xlabel('F','Color',[0 0 1]);
    
    xlabel('F [mL/100mL/min]','Color',[0 0 1]);
    ylabel('E','Color',[0 0 1]);
    zlabel('Vp [mL/100mL]','Color',[0 0 1]);
    title(sprintf('Time resolution = %d sec. Vb = %.2f, Ve = %.2f.',sec_interval,Vb_larss_Murase,Ve_larss_Murase));
    
    %ylabel('Extraction','Color',[0 0 1]);
    
    
    %surf(Real_F_vec,Real_E_vec,Est_Ktrans_vec);
    %scatter(Real_F_vec,Est_Ktrans_vec);
    %scatter3(Real_F_vec,Real_E_vec,Est_Ktrans_vec);
    %ylabel('Estimated Ktrans');
    
    % Print result to PDF
    [Sim_Struct.idx_fig] = Print2Pdf(fig_num, Sim_Struct.idx_fig, 'Vp_vs_E_F.png', './Run_Output/',...
        'Vp vs. E and F', 'VpVsEandF');
    
end


if strcmp(Verbosity,'Full')
    display('-I- Finished Murase/Tofts Comparison...');
end

end