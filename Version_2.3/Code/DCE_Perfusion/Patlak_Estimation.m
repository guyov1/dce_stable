function [ est_Ki_Patlak_noise, est_Vb_Patlak_noise ,E_Patlak_est,  idx_fig ] = Patlak_Estimation( Sim_Struct, est_F_noise, Verbosity, iter_num, avg_num, idx_fig)


% Take from struct variables used in local function
time_vec_minutes                = Sim_Struct.time_vec_minutes;
Sim_AIF_with_noise              = Sim_Struct.Sim_AIF_with_noise;
Sim_Ct_larss_kernel             = Sim_Struct.Sim_Ct_larss_kernel;
plot_flag                       = Sim_Struct.plot_flag;
Ki                              = Sim_Struct.Ki;
Vb_larss                        = Sim_Struct.Vb_larss;
Patlak_Est                      = Sim_Struct.Patlak_Est;

%-----------------------------------
% Ignore points where Ca(t) is very small (because then the result is unstable)
%-----------------------------------
Y_vec_Vb         = Sim_Ct_larss_kernel(:,iter_num,avg_num) ./ Sim_AIF_with_noise(:,iter_num,avg_num); %[mL/100g]
X_vec            = cumtrapz(time_vec_minutes,Sim_AIF_with_noise(:,iter_num,avg_num)) ./ Sim_AIF_with_noise(:,iter_num,avg_num); %[min]
[~, bolus_idx]   = max(diff(Sim_AIF_with_noise(:,iter_num,avg_num)));
% Make sure the index is not too early
if (bolus_idx < 2 )
    base_value       = mean(Sim_AIF_with_noise(1:2,iter_num,avg_num));
else
    base_value       = mean(Sim_AIF_with_noise(1:bolus_idx-1,iter_num,avg_num));
end

mult_val_Thresh  = 3;
Threshold        = mult_val_Thresh*base_value;
stable_idx       = find(Sim_AIF_with_noise(:,iter_num,avg_num) > Threshold);
% Decrease threshold if too many stable points are found
while length(stable_idx) < 2
    mult_val_Thresh  = mult_val_Thresh * 0.8;
    Threshold        = mult_val_Thresh*base_value;
    stable_idx       = find(Sim_AIF_with_noise(:,iter_num,avg_num) > Threshold);
end

% Take the stable points out of the vector
X_vec            = X_vec(stable_idx);
Y_vec_Vb         = Y_vec_Vb(stable_idx);

if (plot_flag)
    
    % plot the AIF and the taken points
    fig_num = figure;
    subplot(2,1,1);
    hold on;
    plot(time_vec_minutes,Sim_AIF_with_noise(:,iter_num,avg_num),'*k');
    plot(time_vec_minutes(stable_idx),Sim_AIF_with_noise(stable_idx,iter_num,avg_num),'og');
    % Plot the line which marks the bolus start
    bolus_start_time = time_vec_minutes(bolus_idx);
    line([bolus_start_time bolus_start_time],[min(Sim_AIF_with_noise(:,iter_num,avg_num)) max(Sim_AIF_with_noise(:,iter_num,avg_num))]);
    plot(time_vec_minutes,Threshold,'-r');
    hold off;
    title('AIF and high value points','fontweight','bold');
    xlabel('[Min]');
    ylabel('AIF [mM]');
    % plot Ct(t) and the taken points
    subplot(2,1,2);
    hold on;
    plot(time_vec_minutes,Sim_Ct_larss_kernel(:,iter_num,avg_num),'*k');
    plot(time_vec_minutes(stable_idx),Sim_Ct_larss_kernel(stable_idx,iter_num,avg_num),'og');
    % Plot the line which marks the bolus start
    bolus_start_time = time_vec_minutes(bolus_idx);
    line([bolus_start_time bolus_start_time],[min(Sim_Ct_larss_kernel(:,iter_num,avg_num)) max(Sim_Ct_larss_kernel(:,iter_num,avg_num))]);
    plot(time_vec_minutes,Threshold,'-r');
    hold off;
    title('Ct(t) and taken points','fontweight','bold');
    xlabel('[Min]');
    ylabel('Ct(t)');
    
    % Print result to PDF
    [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Patlak_Output_1.png', './Run_Output/',...
        'Patlak Results', 'PatlakResults1');
    
end

% Remove Zeros/NaNs/Infs caused because of division by 0
nan_indices           = find(isnan(Y_vec_Vb));
inf_indices           = find(~isfinite(Y_vec_Vb));
Y_vec_Vb(nan_indices) = [];
Y_vec_Vb(inf_indices) = [];
X_vec(nan_indices)    = [];
X_vec(inf_indices)    = [];

% Fine straight line coefficent
[linear_params_spec_points]    = polyfit(X_vec,Y_vec_Vb,1);

% Plot Patlak's fit
a   = linear_params_spec_points(1);
b   = linear_params_spec_points(2);

if (plot_flag)
    
    fig_num = figure;
    subplot(3,1,1);
    hold on;
    vec_size        = max(size(X_vec));
    linear_fit      = linear_params_spec_points(1)*sort(X_vec) + linear_params_spec_points(2);
    %half_idx = round(vec_size/2);
    %scatter(X_vec(1:half_idx),Y_vec(1:half_idx),'g');
    %scatter(X_vec(half_idx+1:end),Y_vec(half_idx+1:end),'r');
    scatter(X_vec,Y_vec_Vb,'g');
    plot(X_vec,linear_fit,'-k');
    h_legend = legend(['Ki=' num2str(Ki,'%.1f') ',Ki-est=' num2str(a,'%.1f')],['Vb=' num2str(Vb_larss,'%.1f') ',Vb-est=' num2str(b,'%.1f')],'Location','NorthWest');
    set(h_legend,'FontSize',7);
    hold off;
    title('Using specified points estimation','fontweight','bold');
    xlabel('\int_{0}^{t} C_a(\tau)d\tau  / C_a(t)');
    ylabel('C_t(t) / C_a(t)');
    %ylabel('$$\int_{0}^{t} \frac{C_t(t)}{C_a(t)}$$', 'interpreter','latex');
    
end

%-----------------------------------
% Using all points result
%-----------------------------------
Y_vec_Vb              = Sim_Ct_larss_kernel(:,iter_num,avg_num) ./ Sim_AIF_with_noise(:,iter_num,avg_num);                           %[mL/100g]
X_vec                 = cumtrapz(time_vec_minutes,Sim_AIF_with_noise(:,iter_num,avg_num)) ./ Sim_AIF_with_noise(:,iter_num,avg_num); %[min]
nan_indices           = find(isnan(Y_vec_Vb));
inf_indices           = find(~isfinite(Y_vec_Vb));
Y_vec_Vb(nan_indices) = [];
Y_vec_Vb(inf_indices) = [];
X_vec(nan_indices)    = [];
X_vec(inf_indices)    = [];

% Fine straight line coefficent
[linear_params_all_points]    = polyfit(X_vec,Y_vec_Vb,1);

% Plot Patlak's fit
a   = linear_params_all_points(1);
b   = linear_params_all_points(2);

if (plot_flag)
    
    subplot(3,1,2);
    hold on;
    vec_size     = max(size(X_vec));
    linear_fit   = linear_params_all_points(1)*sort(X_vec) + linear_params_all_points(2);
    
    scatter(X_vec,Y_vec_Vb,'g');
    plot(sort(X_vec),linear_fit,'-k');
    h_legend = legend(['Ki=' num2str(Ki,'%.1f') ',Ki-est=' num2str(a,'%.1f')],['Vb=' num2str(Vb_larss,'%.1f') ',Vb-est=' num2str(b,'%.1f')],'Location','NorthWest');
    set(h_legend,'FontSize',7);
    hold off;
    title('Using all points estimation','fontweight','bold');
    xlabel('\int_{0}^{t} C_a(\tau)d\tau  / C_a(t)');
    ylabel('C_t(t) / C_a(t)');
    
end

%-----------------------------------
% Linear regression alternative
% Give more weight to the values around the bolus because their SNR is bigger
%-----------------------------------
Y_vec_Vb  = Sim_Ct_larss_kernel(:,iter_num,avg_num) ./ Sim_AIF_with_noise(:,iter_num,avg_num);
X_vec     = cumtrapz(time_vec_minutes,Sim_AIF_with_noise(:,iter_num,avg_num)) ./ Sim_AIF_with_noise(:,iter_num,avg_num);
size_X    = max(size(X_vec));

nan_indices        = find(isnan(Y_vec_Vb));
inf_indices        = find(~isfinite(Y_vec_Vb));
Y_vec_Vb(nan_indices) = [];
Y_vec_Vb(inf_indices) = [];
X_vec(nan_indices) = [];
X_vec(inf_indices) = [];
size_X_new         = max(size(X_vec));

updated_indices    = setdiff(1:size_X,nan_indices);
updated_indices    = setdiff(updated_indices,inf_indices);
relevant_indices   = find(Sim_AIF_with_noise(updated_indices,iter_num,avg_num) > mult_val_Thresh*base_value);

weight_vector                   = ones(1,size_X_new);
weight_value                    = 5;
weight_vector(relevant_indices) = weight_vector(relevant_indices) * weight_value;
W                               = diag(weight_vector);% Weight matrix
A                               = [X_vec ones(size_X_new,1)];
linear_params_weighted_points   = pinv(transpose(A)*W*A)*transpose(A)*W*Y_vec_Vb;

if (plot_flag)
    
    subplot(3,1,3);
    hold on;
    scatter(X_vec,Y_vec_Vb,'k');
    scatter(X_vec(relevant_indices),Y_vec_Vb(relevant_indices),'og');
    linear_fit = linear_params_weighted_points(1)*sort(X_vec) + linear_params_weighted_points(2);
    a          = linear_params_weighted_points(1);
    b          = linear_params_weighted_points(2);
    plot(sort(X_vec),linear_fit,'-k');
    hold off;
    h_legend = legend(['Ki=' num2str(Ki,'%.1f') ',Ki-est=' num2str(a,'%.1f')],['Vb=' num2str(Vb_larss,'%.1f') ',Vb-est=' num2str(b,'%.1f')],'Location','NorthWest');
    set(h_legend,'FontSize',7);
    title('Weighted points','fontweight','bold');
    xlabel('\int_{0}^{t} C_a(\tau)d\tau  / C_a(t)');
    ylabel('C_t(t) / C_a(t)');
    
    % Print result to PDF
    [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Patlak_Output_2.png', './Run_Output/',...
        '', 'PatlakResults2');
    
end

% Choose which Patlak results to take
switch Patlak_Est
    case 'Specified Points'
        est_Ki_Patlak_noise = linear_params_spec_points(1); %a in ax+b
        est_Vb_Patlak_noise = linear_params_spec_points(2); %b in ax+b
    case 'All Points'
        est_Ki_Patlak_noise = linear_params_all_points(1); %a in ax+b
        est_Vb_Patlak_noise = linear_params_all_points(2); %b in ax+b
    case 'Weighted Points'
        est_Ki_Patlak_noise = linear_params_weighted_points(1); %a in ax+b
        est_Vb_Patlak_noise = linear_params_weighted_points(2); %b in ax+b
    otherwise
        error('-E- Patlak estimation not specified correctly!');
end

E_Patlak_est        = est_Ki_Patlak_noise / est_F_noise; % E = Ki / F

end