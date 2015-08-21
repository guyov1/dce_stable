% --------------------- Save old iteration data ---------------------------

Vb_single                   = Sim_Struct.Vb_single;
E_single                    = Sim_Struct.E_single;
F_single                    = Sim_Struct.F_single;
results                     = Sim_Struct.results;

real_larsson_F_vec          = results(15,:);
est_larsson_F_vec           = results(16,:);
error_percent_F             = results(17,:);
std_F                       = results(18,:);

real_t_d_Larss_vec_sec      = results(19,:);
est_t_d_Larss_vec_sec       = results(20,:);
error_percent_t_d_Larss     = results(21,:);
std_t_d_Larss_sec           = results(22,:);

real_larsson_Ktrans_vec          = results(27,:);
est_larsson_Ktrans_vec           = results(28,:);
error_percent_Ktrans             = results(29,:);
std_Ktrans                       = results(30,:);

real_larsson_Vb_vec          = results(39,:);
est_larsson_Vb_vec           = results(40,:);
error_percent_Vb             = results(41,:);
std_Vb                       = results(42,:);


real_larsson_F_vec_last     = real_larsson_F_vec;
real_t_d_Larss_vec_sec_last = real_t_d_Larss_vec_sec;
est_t_d_Larss_vec_sec_last  = est_t_d_Larss_vec_sec;
std_t_d_Larss_sec_last      = std_t_d_Larss_sec;
error_percent_F_last        = error_percent_F;
std_F_last                  = std_F;

real_larsson_Ktrans_vec_last          = real_larsson_Ktrans_vec;
est_larsson_Ktrans_vec_last           = est_larsson_Ktrans_vec;
error_percent_Ktrans_last             = error_percent_Ktrans;
std_Ktrans_last                       = std_Ktrans;

real_larsson_Vb_vec_last          = real_larsson_Vb_vec;
est_larsson_Vb_vec_last           = est_larsson_Vb_vec;
error_percent_Vb_last             = error_percent_Vb;
std_Vb_last                       = std_Vb;

save('./Run_Output/Last_Iter_Vars.mat', 'real_larsson_F_vec_last','real_t_d_Larss_vec_sec_last',...
                                        'est_t_d_Larss_vec_sec_last','std_t_d_Larss_sec_last','error_percent_F_last','std_F_last',...
                                        'real_larsson_Ktrans_vec_last','est_larsson_Ktrans_vec_last','error_percent_Ktrans_last','std_Ktrans_last',...
                                        'real_larsson_Vb_vec_last','est_larsson_Vb_vec_last','error_percent_Vb_last','std_Vb_last');

real_larsson_F_vec_last2     = real_larsson_F_vec;
real_t_d_Larss_vec_sec_last2 = real_t_d_Larss_vec_sec;
est_t_d_Larss_vec_sec_last2  = est_t_d_Larss_vec_sec;
std_t_d_Larss_sec_last2      = std_t_d_Larss_sec;
error_percent_F_last2        = error_percent_F;
std_F_last2                  = std_F;

save('./Run_Output/Last_2_Iter_Vars.mat', 'real_larsson_F_vec_last2', 'real_t_d_Larss_vec_sec_last2','est_t_d_Larss_vec_sec_last2','std_t_d_Larss_sec_last2','error_percent_F_last2','std_F_last2');

real_larsson_F_vec_last3     = real_larsson_F_vec;
real_t_d_Larss_vec_sec_last3 = real_t_d_Larss_vec_sec;
est_t_d_Larss_vec_sec_last3  = est_t_d_Larss_vec_sec;
std_t_d_Larss_sec_last3      = std_t_d_Larss_sec;
error_percent_F_last3        = error_percent_F;
std_F_last3                  = std_F;

save('./Run_Output/Last_3_Iter_Vars.mat', 'real_larsson_F_vec_last3', 'real_t_d_Larss_vec_sec_last3','est_t_d_Larss_vec_sec_last3','std_t_d_Larss_sec_last3','error_percent_F_last3','std_F_last3');

% -------------------------------------------------------------------------

Vb_single                   = Sim_Struct.Vb_single;
E_single                    = Sim_Struct.E_single;
F_single                    = Sim_Struct.F_single;
results                     = Sim_Struct.results;
real_larsson_F_vec          = results(15,:);
est_larsson_F_vec           = results(16,:);
error_percent_F             = results(17,:);
std_F                       = results(18,:);
real_t_d_Larss_vec_sec      = results(19,:);
est_t_d_Larss_vec_sec       = results(20,:);
error_percent_t_d_Larss     = results(21,:);
std_t_d_Larss_sec           = results(22,:);

real_larsson_Ktrans_vec          = results(27,:);
est_larsson_Ktrans_vec           = results(28,:);
error_percent_Ktrans             = results(29,:);
std_Ktrans                       = results(30,:);

real_larsson_Vb_vec          = results(39,:);
est_larsson_Vb_vec           = results(40,:);
error_percent_Vb             = results(41,:);
std_Vb                       = results(42,:);


load('./Run_Output/Last_Iter_Vars.mat');
load('./Run_Output/Last_2_Iter_Vars.mat');
load('./Run_Output/Last_3_Iter_Vars.mat');


%% F error vs AIF delay
fig_num = figure;
%plot(real_sigma_vec,error_percent_sigma);
h1 = errorbar(real_t_d_Larss_vec_sec,error_percent_F,std_F,'b*');

hold on;
h2 = errorbar(real_t_d_Larss_vec_sec_last,error_percent_F_last,std_F_last,'go');
hold off;

title_string = sprintf('est. F error vs. AIF Delay. Vb=%d, E=%.2f, F=%.2f',Vb_single,E_single,real_larsson_F_vec(1));
title(title_string,'FontWeight','bold');
xlabel('AIF delay [Sec]');
ylabel('Error percent');
legend([h1 h2],'No Correction','With Correction'); 

%% Ktrans error vs AIF delay
fig_num = figure;
%plot(real_sigma_vec,error_percent_sigma);
h1 = errorbar(real_t_d_Larss_vec_sec,error_percent_Ktrans,std_Ktrans,'b*');

hold on;
h2 = errorbar(real_t_d_Larss_vec_sec_last,error_percent_Ktrans_last,std_Ktrans_last,'go');
hold off;

title_string = sprintf('est. Ktrans error vs. AIF Delay. Vb=%d, E=%.2f, F=%.2f',Vb_single,E_single,real_larsson_Ktrans_vec(1));
title(title_string,'FontWeight','bold');
xlabel('AIF delay [Sec]');
ylabel('Error percent');
legend([h1 h2],'No Correction','With Correction'); 

%% Vb error vs AIF delay
fig_num = figure;
%plot(real_sigma_vec,error_percent_sigma);
h1 = errorbar(real_t_d_Larss_vec_sec,error_percent_Vb,std_Vb,'b*');

hold on;
h2 = errorbar(real_t_d_Larss_vec_sec_last,error_percent_Vb_last,std_Vb_last,'go');
hold off;

title_string = sprintf('est. Vb error vs. AIF Delay. Vb=%d, E=%.2f, F=%.2f',Vb_single,E_single,real_larsson_Vb_vec(1));
title(title_string,'FontWeight','bold');
xlabel('AIF delay [Sec]');
ylabel('Error percent');
legend([h1 h2],'No Correction','With Correction'); 


%% Estimated vs. real AIF delay

fig_num = figure;

%subplot(2,1,2);
h1 = errorbar(real_t_d_Larss_vec_sec,est_t_d_Larss_vec_sec,std_t_d_Larss_sec,'b');
hr = refline(1); % y=x reference line
set(hr,'Color','k','LineStyle','--','LineWidth',1.5);

hold on;
h2 = errorbar(real_t_d_Larss_vec_sec_last,est_t_d_Larss_vec_sec_last,std_t_d_Larss_sec_last,'g');


% errorbar(real_t_d_Larss_vec_sec_last2,est_t_d_Larss_vec_sec_last2,std_t_d_Larss_sec_last2,'g');
% hr = refline(1); % y=x reference line
% set(hr,'Color','b','LineStyle','--','LineWidth',1.5);

hold off;

title_string = sprintf('Estimated Larss AIF delay vs. original AIF delay. F=%d, Vb=%.2f. (RED-Gaussian)',F_single,Vb_single);
title(title_string,'FontWeight','bold');
xlabel('True AIF delay');
ylabel('Estimated AIF delay');
axis equal;
legend([h1 h2],'Max point time','Novel method'); 


%% -- Real F versus error percent in F estimation

fig_num = figure;

%plot(real_sigma_vec,error_percent_sigma);
h1 = errorbar(real_larsson_F_vec,error_percent_F,std_F,'r');
hold on;
h2 = errorbar(real_larsson_F_vec,error_percent_F_last2,std_F_last2,'b');
h3 = errorbar(real_larsson_F_vec,error_percent_F_last,std_F_last,'g');
hold off;

title_string = sprintf('Est. error vs. original Flow. Vb=%d, E=%.2f, F=%.2f',Vb_single,E_single,real_larsson_F_vec(1));
title(title_string,'FontWeight','bold');
legend([h1 h2 h3],'6 [sec]','4 [sec]','2 [sec]'); 
xlabel('True Flow');
ylabel('Error percent');





