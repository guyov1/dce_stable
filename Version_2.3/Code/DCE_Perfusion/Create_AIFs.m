function [ Sim_Struct ] = Create_AIFs( Sim_Struct, Verbosity )

if strcmp(Verbosity,'Full')
    display('-I- Starting AIFs creation...');
end

% Take from struct variables used in local function
min_interval             = Sim_Struct.min_interval;
High_res_min             = Sim_Struct.High_res_min;
time_vec_minutes         = Sim_Struct.time_vec_minutes;
num_time_points          = length(time_vec_minutes);
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
additional_AIF_delay_min = Sim_Struct.additional_AIF_delay_min;
num_iterations           = Sim_Struct.num_iterations;
num_averages             = Sim_Struct.num_averages;
r_factor                 = Sim_Struct.r_factor;
SNR_ratio                = Sim_Struct.SNR_ratio;


% Initiate  vectors
Sim_Struct.Upsamp_factor            = round( min_interval / High_res_min );
Sim_Struct.Sim_AIF                  = zeros(num_time_points,num_iterations);
Sim_Struct.Sim_AIF_delayed          = zeros(num_time_points,num_iterations);
Sim_Struct.Sim_AIF_high_res         = zeros(num_time_points*Sim_Struct.Upsamp_factor(1),num_iterations);
Sim_Struct.Sim_AIF_delayed_high_res = zeros(num_time_points*Sim_Struct.Upsamp_factor(1),num_iterations);


for iter_num = 1:num_iterations
    
    % Interpolate time vector for highter resolution
    Sim_Struct.time_vec_minutes_high_res = interp(time_vec_minutes,Sim_Struct.Upsamp_factor(iter_num));
    
    % Create the AIF (according to Parker's average population)
    Sim_Struct.Sim_AIF_high_res(:,iter_num)         = AIF_Parker(Sim_Struct.time_vec_minutes_high_res,A1,sig1,T1,A2,sig2,T2,alpha,beta,s,tau,0); %[mM]
    Sim_Struct.Sim_AIF_delayed_high_res(:,iter_num) = AIF_Parker(Sim_Struct.time_vec_minutes_high_res,A1,sig1,T1,A2,sig2,T2,alpha,beta,s,tau,additional_AIF_delay_min(iter_num)); %[mM]
    % Downsample to wanted resolution
    Sim_Struct.Sim_AIF(:,iter_num)                  = downsample(Sim_Struct.Sim_AIF_high_res(:,iter_num),Sim_Struct.Upsamp_factor(iter_num)); %[mM]
    Sim_Struct.Sim_AIF_delayed(:,iter_num)          = downsample(Sim_Struct.Sim_AIF_delayed_high_res(:,iter_num),Sim_Struct.Upsamp_factor(iter_num)); %[mM]
    
    % Normalize to 0->1 and Multiply by relaxivity factor to fit T1 values
    Sim_Struct.Sim_AIF_high_res(:,iter_num)         = r_factor * (Sim_Struct.Sim_AIF_high_res(:,iter_num) / max(Sim_Struct.Sim_AIF_high_res(:,iter_num)) ); %[mM]
    Sim_Struct.Sim_AIF_delayed_high_res(:,iter_num) = r_factor * (Sim_Struct.Sim_AIF_delayed_high_res(:,iter_num) / max(Sim_Struct.Sim_AIF_delayed_high_res(:,iter_num)) ); %[mM]
    Sim_Struct.Sim_AIF(:,iter_num)                  = r_factor * (Sim_Struct.Sim_AIF(:,iter_num) / max(Sim_Struct.Sim_AIF(:,iter_num)) ); %[mM]
    Sim_Struct.Sim_AIF_delayed(:,iter_num)          = r_factor * (Sim_Struct.Sim_AIF_delayed(:,iter_num) / max(Sim_Struct.Sim_AIF_delayed(:,iter_num)) ); %[mM]
    
end

% Adding noise to simulated AIF(t)
Sim_Struct.noise_sigma_AIF_HighRes   = mean(Sim_Struct.Sim_AIF_high_res) ./ SNR_ratio;
Sim_Struct.noise_to_add_AIF_HighRes  = repmat(Sim_Struct.noise_sigma_AIF_HighRes,[num_time_points*Sim_Struct.Upsamp_factor(1) 1 num_averages]) .* randn(num_time_points*Sim_Struct.Upsamp_factor(1),num_iterations,num_averages);

Sim_Struct.noise_sigma_AIF           = mean(Sim_Struct.Sim_AIF) ./ SNR_ratio;
Sim_Struct.noise_to_add_AIF          = repmat(Sim_Struct.noise_sigma_AIF,[num_time_points 1 num_averages]) .* randn(num_time_points,num_iterations,num_averages);

delay_index                          = round(additional_AIF_delay_min ./ min_interval);

% Update the delayed AIF for each iteration
Sim_Struct.noise_to_add_delayed_AIF         = zeros(size(Sim_Struct.noise_to_add_AIF));
Sim_Struct.noise_to_add_delayed_AIF_HighRes = zeros(size(Sim_Struct.noise_to_add_AIF_HighRes));

for iter_num = 1:num_iterations
    for avg_num = 1 : num_averages
        
        % Right shift the AIF
        if (delay_index > 0)
            
            Sim_Struct.noise_to_add_delayed_AIF(1:delay_index(iter_num),iter_num,avg_num)         = Sim_Struct.noise_to_add_AIF(end-delay_index(iter_num)+1:end,iter_num,avg_num);
            Sim_Struct.noise_to_add_delayed_AIF(delay_index(iter_num)+1:end,iter_num,avg_num)     = Sim_Struct.noise_to_add_AIF(1:end-delay_index(iter_num),iter_num,avg_num);
            
            Sim_Struct.noise_to_add_delayed_AIF_HighRes(1:delay_index(iter_num),iter_num,avg_num)         = Sim_Struct.noise_to_add_AIF_HighRes(end-delay_index(iter_num)+1:end,iter_num,avg_num);
            Sim_Struct.noise_to_add_delayed_AIF_HighRes(delay_index(iter_num)+1:end,iter_num,avg_num)     = Sim_Struct.noise_to_add_AIF_HighRes(1:end-delay_index(iter_num),iter_num,avg_num);
            %Sim_Struct.noise_to_add_delayed_AIF   = [Sim_Struct.noise_to_add_AIF(end-delay_index+1:end,:,:) Sim_Struct.noise_to_add_AIF(1:end-delay_index,:,:)];
            
        elseif (delay_index < 0) % Left shift
            
            % Change to a positive index
            delay_index = -1 * delay_index;
            
            Sim_Struct.noise_to_add_AIF(end-delay_index(iter_num)+1:end,iter_num,avg_num) = Sim_Struct.noise_to_add_delayed_AIF(1:delay_index(iter_num),iter_num,avg_num);
            Sim_Struct.noise_to_add_AIF(1:end-delay_index(iter_num),iter_num,avg_num)     = Sim_Struct.noise_to_add_delayed_AIF(delay_index(iter_num)+1:end,iter_num,avg_num);
            
            Sim_Struct.noise_to_add_AIF_HighRes(end-delay_index(iter_num)+1:end,iter_num,avg_num) = Sim_Struct.noise_to_add_delayed_AIF_HighRes(1:delay_index(iter_num),iter_num,avg_num);
            Sim_Struct.noise_to_add_AIF_HighRes(1:end-delay_index(iter_num),iter_num,avg_num)     = Sim_Struct.noise_to_add_delayed_AIF_HighRes(delay_index(iter_num)+1:end,iter_num,avg_num);
            
            %noise_to_add_delayed_AIF   = [Sim_Struct.noise_to_add_AIF(end-delay_index+1:end,:,:) Sim_Struct.noise_to_add_AIF(1:end-delay_index,:,:)];
        else % no shift
            
        end
        
    end
end

% Create AIF matrix
Sim_Struct.Sim_AIF_no_noise           = repmat(Sim_Struct.Sim_AIF, [1 1 num_averages]);
Sim_Struct.Sim_AIF_delayed_no_noise   = repmat(Sim_Struct.Sim_AIF_delayed, [1 1 num_averages]);
Sim_Struct.Sim_AIF_with_noise         = Sim_Struct.Sim_AIF_no_noise         + Sim_Struct.noise_to_add_AIF;
Sim_Struct.Sim_AIF_delayed_with_noise = Sim_Struct.Sim_AIF_delayed_no_noise + Sim_Struct.noise_to_add_delayed_AIF;

Sim_Struct.Sim_AIF_HighRes_no_noise           = repmat(Sim_Struct.Sim_AIF_high_res, [1 1 num_averages]);
Sim_Struct.Sim_AIF_HighRes_delayed_no_noise   = repmat(Sim_Struct.Sim_AIF_delayed_high_res, [1 1 num_averages]);

Sim_Struct.Sim_AIF_HighRes_with_noise         = Sim_Struct.Sim_AIF_HighRes_no_noise         + Sim_Struct.noise_to_add_AIF_HighRes;
Sim_Struct.Sim_AIF_HighRes_delayed_with_noise = Sim_Struct.Sim_AIF_HighRes_delayed_no_noise + Sim_Struct.noise_to_add_delayed_AIF_HighRes;


% Make sure the noise does not create negative values (not physical)
Sim_Struct.Sim_AIF_with_noise         = max(Sim_Struct.Sim_AIF_with_noise,0);
Sim_Struct.Sim_AIF_delayed_with_noise = max(Sim_Struct.Sim_AIF_delayed_with_noise,0);

Sim_Struct.Sim_AIF_HighRes_with_noise         = max(Sim_Struct.Sim_AIF_HighRes_with_noise,0);
Sim_Struct.Sim_AIF_HighRes_delayed_with_noise = max(Sim_Struct.Sim_AIF_HighRes_delayed_with_noise,0);

if strcmp(Verbosity,'Full')
    display('-I- Finished AIFs creation...');
end

end