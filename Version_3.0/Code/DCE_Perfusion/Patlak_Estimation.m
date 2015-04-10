function [ est_Ktrans_Patlak_noise, est_Vb_Patlak_noise ,est_E_Patlak_noise, est_MTT_Patlak_noise, idx_fig ] ...
    = Patlak_Estimation( Sim_Struct, Sim_AIF_with_noise, Sim_Ct_larss_kernel, est_F_noise, Verbosity, iter_num, avg_num, idx_fig)

% Handle zero value
if ( est_F_noise <= 0 )
    ret_value               = NaN;
    est_Ktrans_Patlak_noise = ret_value;
    est_Vb_Patlak_noise     = ret_value;
    est_E_Patlak_noise      = ret_value;
    est_MTT_Patlak_noise    = ret_value;
    return;
end

% Take from struct variables used in local function
time_vec_minutes                = Sim_Struct.time_vec_minutes;
%Sim_AIF_with_noise              = Sim_Struct.Sim_AIF_with_noise;
%Sim_Ct_larss_kernel             = Sim_Struct.Sim_Ct_larss_kernel;
plot_flag                       = Sim_Struct.plot_flag;
Patlak_Est_Type                 = Sim_Struct.Patlak_Est_Type;
Vb_low                          = Sim_Struct.Vb_low;
Ktrans                          = Sim_Struct.Ktrans;       % Simulation ground truth values
Vb_larss                        = Sim_Struct.Vb_larss; % Simulation ground truth values
RealData_Flag                   = Sim_Struct.RealData_Flag;

% Taked needed AIF and CTC from input
if RealData_Flag
    AIF = Sim_AIF_with_noise;
    CTC = Sim_Ct_larss_kernel;
    % Make sure the dimensions agree
    if ~isequal(size(AIF), size(CTC))
        % Transpose one of the vectors
        if size(AIF,1) == 1
            AIF = AIF';
        else
            CTC = CTC';
        end
    end
else
    if sum( size(Sim_AIF_with_noise) ~= 1 ) == 3
        AIF = Sim_AIF_with_noise(:,iter_num,avg_num);
        CTC = Sim_Ct_larss_kernel(:,iter_num,avg_num);
    elseif sum( size(Sim_AIF_with_noise) ~= 1 ) == 2
        AIF = Sim_AIF_with_noise(:,iter_num);
        CTC = Sim_Ct_larss_kernel(:,iter_num);
    else
        AIF = Sim_AIF_with_noise;
        CTC = Sim_Ct_larss_kernel;
    end
end

%-----------------------------------
% Ignore points where Ca(t) is very small (because then the result is unstable)
%-----------------------------------
Y_vec_Vb         = CTC ./ AIF; %[mL/100g]
X_vec            = cumtrapz(time_vec_minutes,AIF) ./ AIF; %[min]
[~, bolus_idx]   = max(diff(AIF));
% Make sure the index is not too early
if (bolus_idx < 2 )
    base_value       = mean(AIF(1:2));
else
    base_value       = mean(AIF(1:bolus_idx-1));
end

mult_val_Thresh  = 3;
Threshold        = mult_val_Thresh*max(0,base_value);
stable_idx       = find(AIF > Threshold);
% Decrease threshold if too many stable points are found
while length(stable_idx) < 2
    mult_val_Thresh  = mult_val_Thresh * 0.8;
    Threshold        = mult_val_Thresh*base_value;
    stable_idx       = find(AIF > Threshold);
end

% Take the stable points out of the vector
X_vec            = X_vec(stable_idx);
Y_vec_Vb         = Y_vec_Vb(stable_idx);

if (plot_flag)
    
    % plot the AIF and the taken points
    fig_num = figure;
    subplot(2,1,1);
    hold on;
    plot(time_vec_minutes,AIF,'*k');
    plot(time_vec_minutes(stable_idx),AIF(stable_idx),'og');
    % Plot the line which marks the bolus start
    bolus_start_time = time_vec_minutes(bolus_idx);
    line([bolus_start_time bolus_start_time],[min(AIF) max(AIF)]);
    plot(time_vec_minutes,Threshold,'-r');
    hold off;
    title('AIF and high value points','fontweight','bold');
    xlabel('[Min]');
    ylabel('AIF [mM]');
    % plot Ct(t) and the taken points
    subplot(2,1,2);
    hold on;
    plot(time_vec_minutes,CTC,'*k');
    plot(time_vec_minutes(stable_idx),CTC(stable_idx),'og');
    % Plot the line which marks the bolus start
    bolus_start_time = time_vec_minutes(bolus_idx);
    line([bolus_start_time bolus_start_time],[min(CTC(:)) max(CTC(:))]);
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

Y_vec_Vb(union(nan_indices,inf_indices)) = [];
X_vec(union(nan_indices,inf_indices))    = [];

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
    h_legend = legend(['Ktrans=' num2str(Ktrans,'%.1f') ',Ktrans-est=' num2str(a,'%.1f')],['Vb=' num2str(Vb_larss,'%.1f') ',Vb-est=' num2str(b,'%.1f')],'Location','NorthWest');
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
Y_vec_Vb              = CTC(:) ./ AIF;                           %[mL/100g]
X_vec                 = cumtrapz(time_vec_minutes,AIF) ./ AIF; %[min]
nan_indices           = find(isnan(Y_vec_Vb));
inf_indices           = find(~isfinite(Y_vec_Vb));
Y_vec_Vb(union(nan_indices,inf_indices)) = [];
X_vec(union(nan_indices,inf_indices))    = [];

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
    h_legend = legend(['Ktrans=' num2str(Ktrans,'%.1f') ',Ktrans-est=' num2str(a,'%.1f')],['Vb=' num2str(Vb_larss,'%.1f') ',Vb-est=' num2str(b,'%.1f')],'Location','NorthWest');
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
Y_vec_Vb  = CTC(:) ./ AIF;
X_vec     = cumtrapz(time_vec_minutes,AIF) ./ AIF;
size_X    = max(size(X_vec));

nan_indices                              = find(isnan(Y_vec_Vb));
inf_indices                              = find(~isfinite(Y_vec_Vb));
Y_vec_Vb(union(nan_indices,inf_indices)) = [];
X_vec(union(nan_indices,inf_indices))    = [];
size_X_new                               = max(size(X_vec));

updated_indices    = setdiff(1:size_X,nan_indices);
updated_indices    = setdiff(updated_indices,inf_indices);
relevant_indices   = find(AIF(updated_indices) > mult_val_Thresh*base_value);

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
    h_legend = legend(['Ktrans=' num2str(Ktrans,'%.1f') ',Ktrans-est=' num2str(a,'%.1f')],['Vb=' num2str(Vb_larss,'%.1f') ',Vb-est=' num2str(b,'%.1f')],'Location','NorthWest');
    set(h_legend,'FontSize',7);
    title('Weighted points','fontweight','bold');
    xlabel('\int_{0}^{t} C_a(\tau)d\tau  / C_a(t)');
    ylabel('C_t(t) / C_a(t)');
    
    % Print result to PDF
    [idx_fig] = Print2Pdf(fig_num, idx_fig, 'Patlak_Output_2.png', './Run_Output/',...
        '', 'PatlakResults2');
    
end

% Choose which Patlak results to take
switch Patlak_Est_Type
    case 'Specified Points'
        est_Ktrans_Patlak_noise = max(0, linear_params_spec_points(1)); %a in ax+b
        est_Vb_Patlak_noise = max(Vb_low, linear_params_spec_points(2)); %b in ax+b
        
        %est_Ktrans_Patlak_noise = linear_params_spec_points(1); %a in ax+b
        %est_Vb_Patlak_noise = linear_params_spec_points(2); %b in ax+b
    case 'All Points'
        est_Ktrans_Patlak_noise = max(0, linear_params_all_points(1)); %a in ax+b
        est_Vb_Patlak_noise = max(Vb_low, linear_params_all_points(2)); %b in ax+b
    case 'Weighted Points'
        est_Ktrans_Patlak_noise = max(0, linear_params_weighted_points(1)); %a in ax+b
        est_Vb_Patlak_noise = max(Vb_low, linear_params_weighted_points(2)); %b in ax+b
    otherwise
        error('-E- Patlak estimation not specified correctly!');
end

est_E_Patlak_noise    = est_Ktrans_Patlak_noise / est_F_noise; % E = Ktrans / F
est_MTT_Patlak_noise  = est_Vb_Patlak_noise / est_F_noise;

end