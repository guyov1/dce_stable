
% Take data saved from different simulation
for sim_iter = 1:4
    eval(['Data = load(''./For_Article_2/Saved_iter_' num2str(sim_iter) '.mat'')']);
    
    eval(['results' num2str(sim_iter) '= Data;']);

end

for sim_iter = 1:4
        
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
   
end

Vb_single  = 12;
E_single   = 0.1;

% Vb_single = Vb_single.Vb_single;
% E_single  = E_single.E_single;

fig_num = figure;

subplot(2,1,1);
%plot(real_sigma_vec,error_percent_sigma);

h1 = errorbar(real_larsson_F_vec1,error_percent_F1,std_F1,'k');

hold on;
h2 = errorbar(real_larsson_F_vec2,error_percent_F2,std_F2,'b');
h3 = errorbar(real_larsson_F_vec3,error_percent_F3,std_F3,'g');
h4 = errorbar(real_larsson_F_vec4,error_percent_F4,std_F4,'m');
%h5 = errorbar(real_larsson_F_vec1,error_percent_Sourbron_F1,std_Sourbron_F1,'r');
hold off;



title_string = sprintf('est. error vs. original Flow. Vb=%d, E=%.2f',Vb_single,E_single);
title(title_string,'FontWeight','bold');
xlabel('Original Flow');
ylabel('Error percent');
legend([h1 h2 h3 h4],'PCA_2nd','Splines 1st deriv','Splines 2nd deriv','Wiener'); 

subplot(2,1,2);

h1 = errorbar(real_larsson_F_vec1,est_larsson_F_vec1,std_F1,'k');
hold on;
h2 = errorbar(real_larsson_F_vec2,est_larsson_F_vec2,std_F2,'b');
h3 = errorbar(real_larsson_F_vec3,est_larsson_F_vec3,std_F3,'g');
h4 = errorbar(real_larsson_F_vec4,est_larsson_F_vec4,std_F4,'m');
%h5 = errorbar(real_larsson_F_vec1,est_larsson_F_Sourbron_vec1,std_Sourbron_F1,'r');

hold off;

hr = refline(1); % y=x reference line
set(hr,'Color','k','LineStyle','--','LineWidth',1.5);

title_string = sprintf('est. Flow vs. original Flow. Vb=%d, E=%.2f',Vb_single,E_single);
title(title_string,'FontWeight','bold');
xlabel('Original Flow');
ylabel('est. Flow');
legend([h1 h2 h3 h4],'PCA','Splines 1st deriv','Splines 2nd deriv','Wiener'); 

