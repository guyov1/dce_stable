function [ Sim_Struct ] = Driving_Differential_Eq( Sim_Struct, Verbosity )


if ~strcmp(Verbosity,'None')
    display('-I- Starting driving differential equation...');
end

% Take from struct variables used in local function
min_interval             = Sim_Struct.min_interval;
time_vec_minutes         = Sim_Struct.time_vec_minutes;
A1                       = Sim_Struct.A1;
sig1                     = Sim_Struct.sig1;
T1                       = Sim_Struct.T1;
A2                       = Sim_Struct.A2;
sig2                     = Sim_Struct.sig2;
T2                       = Sim_Struct.T2;
alpha                    = Sim_Struct.alpha;
beta                     = Sim_Struct.beta;
s                        = Sim_Struct.s;
tau                      = Sim_Struct.tau;
r_factor                 = Sim_Struct.r_factor;
total_sim_time_min       = Sim_Struct.total_sim_time_min;
Hct                      = Sim_Struct.Hct;
F                        = Sim_Struct.F;
Ktrans                   = Sim_Struct.Ktrans;
Ve_larss                 = Sim_Struct.Ve_larss;
Vb_larss                 = Sim_Struct.Vb_larss;
E                        = Sim_Struct.E;
adjusted_larsson         = Sim_Struct.Adjusted_Larsson_Model;

% Iteration index to take for calculation
iter_idx = 1;

% Initiating high resolution time vectors
sec_interval_HighRes       = 0.01;                %[sec]
min_interval_HighRes       = sec_interval_HighRes ./ 60;  %[min]
num_time_stamps_test       = round(total_sim_time_min ./ min_interval_HighRes);
time_vec_minutes_HighRes   = (0:num_time_stamps_test-1).* min_interval_HighRes;

% Creating high resolution AIF
Sim_AIF_HighRes            = AIF_Parker(time_vec_minutes_HighRes,A1,sig1,T1,A2,sig2,T2,alpha,beta,s,tau,0); %[mM]
Sim_AIF_HighRes            = Sim_AIF_HighRes / max(Sim_AIF_HighRes);
Sim_AIF_HighRes            = Sim_AIF_HighRes * r_factor;

% Driving the differential equation to get Cb(nT) and Ce(nT)
[ Cb_vec, Ce_vec ]         = Larsson_Differential_Equation( time_vec_minutes_HighRes, Vb_larss(iter_idx),...
                                                            Ve_larss(iter_idx), Sim_AIF_HighRes, Hct(iter_idx),...
                                                            F(iter_idx), Ktrans(iter_idx));
% Creating high resolution larsson filter
if (adjusted_larsson)
    IRF_larss_HighRes          = Adjusted_Larsson_Filter(time_vec_minutes_HighRes, F(iter_idx),...
        Vb_larss(iter_idx), E(iter_idx), Ve_larss(iter_idx));  % No units
else
    IRF_larss_HighRes          = Larsson_Filter(time_vec_minutes_HighRes, F(iter_idx),...
        Vb_larss(iter_idx), E(iter_idx), Ve_larss(iter_idx), Hct(iter_idx));  % No units
end

larss_filter_HighRes       = F(iter_idx)*IRF_larss_HighRes; % [mL/100g/min]
% Computing C(t) using the regular convolution with double exponentials
if Sim_Struct.ignore_time_delta
    Sim_Ct_integral_res        = filter(Sim_Struct.larss_filter(:,iter_idx),1,Sim_Struct.Sim_AIF_delayed_no_noise);
else
    Sim_Ct_integral_res        = filter(Sim_Struct.larss_filter(:,iter_idx)*min_interval(iter_idx),1,Sim_Struct.Sim_AIF_delayed_no_noise);
end

% Computing Ct(nT) out of Cb(nT) and Ce(nT)
Ct_vec                     = Vb_larss(iter_idx)*Cb_vec + Ve_larss(iter_idx)*Ce_vec;

if (Sim_Struct.plot_flag)
    
    fig_num = figure;
    subplot(4,1,1);
    hold on;
    h1 = plot(time_vec_minutes_HighRes,Sim_AIF_HighRes,'*g');
    h2 = plot(time_vec_minutes,Sim_Struct.Sim_AIF_no_noise(:,iter_idx,iter_idx),'ob');
    h3 = plot(time_vec_minutes_HighRes,Ct_vec,'*m');
    h4 = plot(time_vec_minutes,Sim_Ct_integral_res(:,iter_idx,iter_idx),'or');
    
    hold off;
    title('AIF and Ct','FontWeight','bold');
    xlabel('Time [Min]');
    legend([h1 h2 h3 h4],'Differential Driven AIF','Simulated AIF','Differential Driven Ct','Simulated Ct');
    
    subplot(4,1,2);
    hold on;
    h1 = plot(time_vec_minutes_HighRes,Cb_vec,'*');
    h2 = plot(time_vec_minutes_HighRes,Ce_vec,'g*');
    hold off;
    title('Differential driven Cb(t) and Ce(t)','FontWeight','bold');
    xlabel('Time [Min]');
    legend([h1 h2],'Differential Driven Cb','Differential Driven Ce');
    
    subplot(4,1,3);
    hold on;
    h1 = plot(time_vec_minutes_HighRes,Vb_larss(iter_idx)*Cb_vec,'*');
    h2 = plot(time_vec_minutes_HighRes,Ve_larss(iter_idx)*Ce_vec,'g*');
    hold off;
    title('Differential driven Vb*Cb(t) and Ve*Ce(t)','FontWeight','bold');
    xlabel('Time [Min]');
    legend([h1 h2],'Differential Driven Vb*Cb','Differential Driven Ve*Ce');
    
    subplot(4,1,4);
    hold on;
    h1 = plot(time_vec_minutes_HighRes,Ct_vec,'*b');
    h2 = plot(time_vec_minutes,Sim_Ct_integral_res(:,iter_idx,iter_idx),'og');
    hold off;
    title('Simulated Ct Vs. Differential driven Ct','FontWeight','bold');
    xlabel('Time [Min]');
    legend([h1 h2],'Differential Driven Ct','Simulated Ct');
    
    % Print result to PDF
    [Sim_Struct.idx_fig] = Print2Pdf(fig_num, Sim_Struct.idx_fig, 'Sim_Ct_Vs_Differential_Result.png', './Run_Output/',...
        'Sim Ct Vs. Differential driven', 'DifferentialCt');
end

if strcmp(Verbosity,'Full')
    display('-I- Finished driving differential equation...');
end

end