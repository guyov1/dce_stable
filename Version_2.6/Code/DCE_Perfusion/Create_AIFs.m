function [ Sim_Struct ] = Create_AIFs( Sim_Struct, Verbosity )

tic; 

if ~strcmp(Verbosity,'None')
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
AIF_delay_low            = Sim_Struct.AIF_delay_low;
AIF_delay_max            = Sim_Struct.AIF_delay_max;
num_iterations           = Sim_Struct.num_iterations;
num_averages             = Sim_Struct.num_averages;
r_factor                 = Sim_Struct.r_factor;
SNR_ratio                = Sim_Struct.SNR_ratio;

if (Sim_Struct.Add_Randomly_AIF_Delay)
    Sim_Struct.additional_AIF_delay_sec = AIF_delay_low  + (AIF_delay_max  - AIF_delay_low)  * rand(1, num_iterations);
    Sim_Struct.additional_AIF_delay_min = Sim_Struct.additional_AIF_delay_sec / 60;
    additional_AIF_delay_min            = Sim_Struct.additional_AIF_delay_min;
end

% Initiate  vectors
Sim_Struct.Upsamp_factor = round( min_interval / High_res_min );
Sim_AIF                  = zeros(num_time_points,num_iterations);
Sim_AIF_delayed          = zeros(num_time_points,num_iterations);

% CURRENTLY NOT SUPPORYING MULTIPLE UPSAMPLING FACTORS
Upsamp_factor = Sim_Struct.Upsamp_factor(1);

Sim_AIF_high_res         = zeros(num_time_points*Upsamp_factor,num_iterations);
Sim_AIF_delayed_high_res = zeros(num_time_points*Upsamp_factor,num_iterations);

time_vec_minutes_high_res = interp( time_vec_minutes, Upsamp_factor ); % Interpolate time vector for highter resolution

%Sim_Struct.time_vec_minutes_high_res = interp(time_vec_minutes,Sim_Struct.Upsamp_factor(iter_num)); % Interpolate time vector for highter resolution

if ~Sim_Struct.FORCE_SERIAL   
    parfor iter_num = 1:num_iterations
        [ Sim_AIF_high_res(:,iter_num), Sim_AIF_delayed_high_res(:,iter_num), Sim_AIF(:,iter_num), Sim_AIF_delayed(:,iter_num)] ...
            = Create_Wanted_AIFs(additional_AIF_delay_min(iter_num), time_vec_minutes_high_res,A1,sig1,T1,A2,sig2,T2,alpha,beta,s,tau, r_factor, Upsamp_factor );
        
        if ( mod(iter_num,100) == 0 )
            display(sprintf('-I- Finished creating 100 AIFs...'));
        end
        
    end
else
    for iter_num = 1:num_iterations
        [ Sim_AIF_high_res(:,iter_num), Sim_AIF_delayed_high_res(:,iter_num), Sim_AIF(:,iter_num), Sim_AIF_delayed(:,iter_num)] ...
            = Create_Wanted_AIFs(additional_AIF_delay_min(iter_num), time_vec_minutes_high_res,A1,sig1,T1,A2,sig2,T2,alpha,beta,s,tau, r_factor, Upsamp_factor );
        
        if ( mod(iter_num,100) == 0 )
            display(sprintf('-I- Finished creating %d AIFs...',iter_num));
        end
        
    end
end

% Adding noise to simulated AIF(t)
Sim_Struct.noise_sigma_AIF_HighRes   = mean(Sim_AIF_high_res) ./ SNR_ratio;
Sim_Struct.noise_to_add_AIF_HighRes  = repmat(Sim_Struct.noise_sigma_AIF_HighRes,[num_time_points*Upsamp_factor 1 num_averages]) .* randn(num_time_points*Upsamp_factor,num_iterations,num_averages);

Sim_Struct.noise_sigma_AIF           = mean(Sim_AIF) ./ SNR_ratio;
Sim_Struct.noise_to_add_AIF          = repmat(Sim_Struct.noise_sigma_AIF,[num_time_points 1 num_averages]) .* randn(num_time_points,num_iterations,num_averages);

delay_index                          = round(additional_AIF_delay_min ./ min_interval);

% Update the delayed AIF noise for each iteration
[ Sim_Struct.noise_to_add_delayed_AIF, Sim_Struct.noise_to_add_delayed_AIF_HighRes] = ...
    Update_Delayed_AIF_Noise(Sim_Struct.noise_to_add_AIF, Sim_Struct.noise_to_add_AIF_HighRes, delay_index, Sim_Struct.FORCE_SERIAL, num_iterations, num_averages);

% Create AIF matrix
Sim_Struct.Sim_AIF_no_noise           = repmat(Sim_AIF, [1 1 num_averages]);
Sim_Struct.Sim_AIF_delayed_no_noise   = repmat(Sim_AIF_delayed, [1 1 num_averages]);
Sim_Struct.Sim_AIF_with_noise         = Sim_Struct.Sim_AIF_no_noise         + Sim_Struct.noise_to_add_AIF;
Sim_Struct.Sim_AIF_delayed_with_noise = Sim_Struct.Sim_AIF_delayed_no_noise + Sim_Struct.noise_to_add_delayed_AIF;

Sim_Struct.Sim_AIF_HighRes_no_noise           = repmat(Sim_AIF_high_res, [1 1 num_averages]);
Sim_Struct.Sim_AIF_HighRes_delayed_no_noise   = repmat(Sim_AIF_delayed_high_res, [1 1 num_averages]);

Sim_Struct.Sim_AIF_HighRes_with_noise         = Sim_Struct.Sim_AIF_HighRes_no_noise         + Sim_Struct.noise_to_add_AIF_HighRes;
Sim_Struct.Sim_AIF_HighRes_delayed_with_noise = Sim_Struct.Sim_AIF_HighRes_delayed_no_noise + Sim_Struct.noise_to_add_delayed_AIF_HighRes;


% Make sure the noise does not create negative values (not physical)
Sim_Struct.Sim_AIF_with_noise         = max(Sim_Struct.Sim_AIF_with_noise,0);
Sim_Struct.Sim_AIF_delayed_with_noise = max(Sim_Struct.Sim_AIF_delayed_with_noise,0);

Sim_Struct.Sim_AIF_HighRes_with_noise         = max(Sim_Struct.Sim_AIF_HighRes_with_noise,0);
Sim_Struct.Sim_AIF_HighRes_delayed_with_noise = max(Sim_Struct.Sim_AIF_HighRes_delayed_with_noise,0);

% Update output struct
Sim_Struct.time_vec_minutes_high_res = time_vec_minutes_high_res;
Sim_Struct.Sim_AIF_high_res          = Sim_AIF_high_res;
Sim_Struct.Sim_AIF_delayed_high_res  = Sim_AIF_delayed_high_res;
Sim_Struct.Sim_AIF                   = Sim_AIF;
Sim_Struct.Sim_AIF_delayed           = Sim_AIF_delayed;

time_finish = toc;
if ~strcmp(Verbosity,'None')
    display(sprintf('-I- Creating AIFs took %.2f seconds to finish...',time_finish));
end


if strcmp(Verbosity,'Full')
    display('-I- Finished AIFs creation...');
end

end