
% Take data saved from different simulation
for sim_iter = 1:2
    eval(['Data = load(''./For_Article_2/Saved_iter_' num2str(sim_iter) '.mat'')']);
    
    eval(['results' num2str(sim_iter) '= Data;']);
    
end

for sim_iter = 1:2
    
    % F
    eval(['real_larsson_F_vec' num2str(sim_iter) '= results' num2str(sim_iter) '.results(15,:);']);
    eval(['est_larsson_F_vec' num2str(sim_iter) '= results' num2str(sim_iter) '.results(16,:);']);
    eval(['error_percent_F' num2str(sim_iter) '= results' num2str(sim_iter) '.results(17,:);']);
    eval(['std_F' num2str(sim_iter) '= results' num2str(sim_iter) '.results(18,:);']);
    
    
    % Sourbron - F
    eval(['real_larsson_F_Sourbron_vec' num2str(sim_iter) '= results' num2str(sim_iter) '.results(67,:);']);
    eval(['est_larsson_F_Sourbron_vec' num2str(sim_iter) '= results' num2str(sim_iter) '.results(68,:);']);
    eval(['error_percent_Sourbron_F' num2str(sim_iter) '= results' num2str(sim_iter) '.results(69,:);']);
    eval(['std_Sourbron_F' num2str(sim_iter) '= results' num2str(sim_iter) '.results(70,:);']);
    
    
    eval(['real_t_d_Larss_vec_sec' num2str(sim_iter) '= results' num2str(sim_iter) '.results(19,:);']);
    eval(['est_t_d_Larss_vec_sec' num2str(sim_iter) '= results' num2str(sim_iter) '.results(20,:);']);
    eval(['error_percent_t_d_Larss' num2str(sim_iter) '= results' num2str(sim_iter) '.results(21,:);']);
    eval(['std_t_d_Larss_sec' num2str(sim_iter) '= results' num2str(sim_iter) '.results(22,:);']);
    
    
end

%Vb_single  = 12;
%E_single   = 0.1;
F_single   = 60;

Vb_single = results1.Vb_single;
E_single  = results1.E_single;
%F_single  = results1.F_single;

fig_num = figure;

%subplot(2,1,2);
h1 = errorbar(real_t_d_Larss_vec_sec1,est_t_d_Larss_vec_sec1,std_t_d_Larss_sec1,'b');
hold on;
h2 = errorbar(real_t_d_Larss_vec_sec2,est_t_d_Larss_vec_sec2,std_t_d_Larss_sec2,'g');
hold off;
hr = refline(1); % y=x reference line
set(hr,'Color','k','LineStyle','--','LineWidth',1.5);

title_string = sprintf('est. Larss AIF delay vs. original AIF delay. F=%d, Vb=%.2f.',F_single,Vb_single);
title(title_string,'FontWeight','bold');
xlabel('True AIF delay');
ylabel('est. AIF delay');
legend([h1 h2],'Max IRF','BiExp method');


% ---------------------------------------------------------------------------------------------

fig_num = figure;

subplot(2,1,1);
%plot(real_sigma_vec,error_percent_sigma);

h1 = errorbar(real_larsson_F_vec1,error_percent_F1,std_F1,'k');

hold on;
h2 = errorbar(real_larsson_F_vec2,error_percent_F2,std_F2,'g');
%h5 = errorbar(real_larsson_F_vec1,error_percent_Sourbron_F1,std_Sourbron_F1,'r');
hold off;



title_string = sprintf('est. error vs. original Flow. Vb=%d, E=%.2f',Vb_single,E_single);
title(title_string,'FontWeight','bold');
xlabel('True Flow');
ylabel('Error percent');
legend([h1 h2],'No correction','With correction'); 

subplot(2,1,2);

h1 = errorbar(real_larsson_F_vec1,est_larsson_F_vec1,std_F1,'k');
hold on;
h2 = errorbar(real_larsson_F_vec2,est_larsson_F_vec2,std_F2,'g');

%h5 = errorbar(real_larsson_F_vec1,est_larsson_F_Sourbron_vec1,std_Sourbron_F1,'r');

hold off;

hr = refline(1); % y=x reference line
set(hr,'Color','k','LineStyle','--','LineWidth',1.5);

title_string = sprintf('est. Flow vs. original Flow. Vb=%d, E=%.2f',Vb_single,E_single);
title(title_string,'FontWeight','bold');
xlabel('True Flow');
ylabel('est. Flow');
legend([h1 h2],'No correction','With correction'); 