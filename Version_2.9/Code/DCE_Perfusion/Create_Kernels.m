function [ Sim_Struct ] = Create_Kernels( Sim_Struct, Verbosity )

% ------------------------------------------------------------------
% ------------ Create_Kernels --------------------------------------
% ------------------------------------------------------------------
% ------ Creates Gaussian, Larsson, Sourbron and Patlak IRFs -------
% ------ duplicated to number of iterations                  -------
% ------------------------------------------------------------------


tic;

if ~strcmp(Verbosity,'None')
    display('-I- Starting Kernels creation...');
end

% Take from struct variables used in local function
time_vec_minutes          = Sim_Struct.time_vec_minutes;
num_iterations            = Sim_Struct.num_iterations;
time_vec_minutes_high_res = Sim_Struct.time_vec_minutes_high_res;
adjusted_larsson          = Sim_Struct.Adjusted_Larsson_Model;
F_low                     = Sim_Struct.F_low;
F_max                     = Sim_Struct.F_max;
Vb_low                    = Sim_Struct.Vb_low;
Vb_max                    = Sim_Struct.Vb_max;
Ve_low                    = Sim_Struct.Ve_low;
Ve_max                    = Sim_Struct.Ve_max;
E_low                     = Sim_Struct.E_low;
E_max                     = Sim_Struct.E_max;

% ------------------------
% Gaussian filter
% ------------------------
Sim_Struct.gauss_filter         = zeros(size(time_vec_minutes,2),num_iterations);
Sim_Struct.gauss_filter_HighRes = zeros(size(time_vec_minutes_high_res,2),num_iterations);

for i=1:num_iterations
    Sim_Struct.gauss_filter(:,i)         = Gaussian(time_vec_minutes, Sim_Struct.t_d(i), Sim_Struct.sigma(i)^2, Sim_Struct.amplitude(i));
    Sim_Struct.gauss_filter_HighRes(:,i) = Gaussian(time_vec_minutes_high_res, Sim_Struct.t_d(i), Sim_Struct.sigma(i)^2, Sim_Struct.amplitude(i));
end

% ------------------------
% Larrson Filter
% as described in Larsson's article
% ------------------------

if (Sim_Struct.iterate_F_larsson)
    Sim_Struct.F = Sim_Struct.F_vec;
else
    Sim_Struct.F = repmat(Sim_Struct.F_single,1,num_iterations); % [mL/100g/min] , SNR -> 7. If F=60, SNR->17
end

if (Sim_Struct.iterate_Vb_larsson)
    Sim_Struct.Vb_larss =  Sim_Struct.Vb_vec;
else
    Sim_Struct.Vb_larss = repmat(Sim_Struct.Vb_single,1,num_iterations);       % [mL/100g]     , they used 3,6,12,18
end

if (Sim_Struct.iterate_E_larsson)
    Sim_Struct.E = Sim_Struct.E_vec;
else
    Sim_Struct.E = repmat(Sim_Struct.E_single,1,num_iterations);     % Extraction. They used 0,0.1,0.5,1
end

if (Sim_Struct.iterate_Ve_larsson)
    Sim_Struct.Ve_larss = Sim_Struct.Ve_vec;
else
    Sim_Struct.Ve_larss = Sim_Struct.Ve_single * ( 1 - (Sim_Struct.Vb_larss/100) )*100; % Must be smaller than Vtis
end

if (Sim_Struct.iterate_uniformly)
    Sim_Struct.F        = F_low  + (F_max  - F_low)  * rand(1, num_iterations);
    Sim_Struct.E        = E_low  + (E_max  - E_low)  * rand(1, num_iterations);
    Sim_Struct.Vb_larss = Vb_low + (Vb_max - Vb_low) * rand(1, num_iterations);
    Sim_Struct.Ve_larss = Ve_low + (Ve_max - Ve_low) * rand(1, num_iterations);
end


% Hematocrit according to Larsson's article
Sim_Struct.Hct          = repmat(Sim_Struct.Hct_single,1,num_iterations);

Sim_Struct.IRF_larss         = zeros(size(time_vec_minutes,2),num_iterations);
Sim_Struct.IRF_larss_HighRes = zeros(size(time_vec_minutes_high_res,2),num_iterations);
for i=1:num_iterations
    
    if (adjusted_larsson)
        Sim_Struct.IRF_larss(:,i)         = Adjusted_Larsson_Filter(time_vec_minutes, Sim_Struct.F(i), Sim_Struct.Vb_larss(i), Sim_Struct.E(i),...
            Sim_Struct.Ve_larss(i));  % No units
        Sim_Struct.IRF_larss_HighRes(:,i) = Adjusted_Larsson_Filter(time_vec_minutes_high_res, Sim_Struct.F(i), Sim_Struct.Vb_larss(i), Sim_Struct.E(i),...
            Sim_Struct.Ve_larss(i));  % No units
        
    else
        Sim_Struct.IRF_larss(:,i)         = Larsson_Filter(time_vec_minutes, Sim_Struct.F(i), Sim_Struct.Vb_larss(i), Sim_Struct.E(i),...
            Sim_Struct.Ve_larss(i), Sim_Struct.Hct(i));  % No units
        Sim_Struct.IRF_larss_HighRes(:,i) = Larsson_Filter(time_vec_minutes_high_res, Sim_Struct.F(i), Sim_Struct.Vb_larss(i), Sim_Struct.E(i),...
            Sim_Struct.Ve_larss(i), Sim_Struct.Hct(i));  % No units
    end
    
    
end


Sim_Struct.larss_filter           = repmat(Sim_Struct.F,[size(time_vec_minutes,2) 1]) .* Sim_Struct.IRF_larss; % [mL/100g/min]
Sim_Struct.larss_filter_HighRes   = repmat(Sim_Struct.F,[size(time_vec_minutes_high_res,2) 1]) .* Sim_Struct.IRF_larss_HighRes; % [mL/100g/min]
Sim_Struct.CBF                    = Sim_Struct.F;
Sim_Struct.CBV                    = Sim_Struct.Vb_larss;
Sim_Struct.Vd                     = Sim_Struct.Vb_larss + Sim_Struct.Ve_larss;

% MTT equals integration over IRF
Sim_Struct.MTT                    = cumtrapz(time_vec_minutes,Sim_Struct.IRF_larss);
Sim_Struct.MTT                    = Sim_Struct.MTT(end,:);
Sim_Struct.Vd                     = Sim_Struct.Vb_larss + Sim_Struct.Ve_larss; %Vd = F * MTT;
Sim_Struct.Ktrans                     = Sim_Struct.E .* Sim_Struct.F;

%Sim_Struct.PS                     = -Sim_Struct.F .* log(1-(Sim_Struct.Ktrans ./ Sim_Struct.F));
Sim_Struct.PS                     = ( Sim_Struct.E .* Sim_Struct.F ) ./ (1 - Sim_Struct.E);                             % [mL/100g/min]

% ------------------------
% Sourbron Filter
% as described in Sourbron's article
% ------------------------
Sim_Struct.Fp                = Sim_Struct.F;                            % [mL/100g/min], Same as Fp in article
Sim_Struct.Fe                = Sim_Struct.PS; %Sim_Struct.Ktrans
Sim_Struct.Vp_sourbron       = Sim_Struct.Vb_larss;
Sim_Struct.Ve_sourbron       = Sim_Struct.Ve_larss; % Must be smaller than Vtis

Sim_Struct.IRF_sourbron = zeros(size(time_vec_minutes,2),num_iterations);
for i=1:num_iterations
    Sim_Struct.IRF_sourbron(:,i)      = Sourbron_Filter(time_vec_minutes, Sim_Struct.Fp(i), Sim_Struct.Vp_sourbron(i),...
                                                        Sim_Struct.Fe(i), Sim_Struct.Ve_sourbron(i));  % No units
end

Sim_Struct.sourbron_filter   = repmat(Sim_Struct.Fp,[size(time_vec_minutes,2) 1]) .* Sim_Struct.IRF_sourbron; % [mL/100g/min]
%larss_filter      = sourbron_filter;

% Plotting sourbron vs. larsson filter
if (Sim_Struct.plot_flag)
    
    fig_num = figure;
    hold on;
    h1 = plot(time_vec_minutes,Sim_Struct.sourbron_filter,'*r');
    h2 = plot(time_vec_minutes,Sim_Struct.larss_filter,'og');
    hold off;
    title('Larsson Vs. Sourbron filter','FontWeight','bold');
    xlabel('Time [Min]');
    legend([h1(1) h2(1)],'Sourbron Filter','Larsson Filter');
    
    % Display F values
    annotation('textbox',...
        [0.6 0.49 0.1 0.1],...
        'String',{['F = ' num2str(Sim_Struct.F) '  [mL/100g]']},...
        'FontSize',8,...
        'FontName','Arial',...
        'LineStyle','-',...
        'EdgeColor',[0 0 0],...
        'LineWidth',2,...
        'BackgroundColor',[1 1 1],...
        'Color',[0.0 0.0 0]);
    
    % Print result to PDF
    [Sim_Struct.idx_fig] = Print2Pdf(fig_num, Sim_Struct.idx_fig, 'Larsson_vs_Sourbron.png', './Run_Output/', 'Larsson Vs. Sourbron filter', 'LarssonVsSourbron');
    
end

% ------------------------
% Patlak Filter
% ------------------------
Sim_Struct.vb_delta               = zeros(size(Sim_Struct.larss_filter));
Sim_Struct.vb_delta(1,:)          = Sim_Struct.Vb_larss;
Sim_Struct.Patlak_filter          = Sim_Struct.vb_delta + repmat(Sim_Struct.Ktrans, [size(Sim_Struct.larss_filter,1) 1]); % [mL/100g/min]

time_finish = toc;
if ~strcmp(Verbosity,'None')
    display(sprintf('-I- Creating Kernels took %.2f seconds to finish...',time_finish));
end

if strcmp(Verbosity,'Full')
    display('-I- Finished kernels creation...');
end

end