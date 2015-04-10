function [Sim_Ct_gauss_kernel, Sim_Ct_larss_kernel, noise_to_add_gauss, noise_to_add_larss]= Filter_AIF_Process(Sim_Struct, Sim_Ct_gauss_kernel, Sim_Ct_larss_kernel, noise_to_add_gauss, noise_to_add_larss);

% Take from struct variables used in local function
min_interval                     = Sim_Struct.min_interval;

Sim_AIF_HighRes_delayed_no_noise = Sim_Struct.Sim_AIF_HighRes_delayed_no_noise;
num_iterations                   = Sim_Struct.num_iterations;
num_averages                     = Sim_Struct.num_averages;
SNR_ratio                        = Sim_Struct.SNR_ratio;
High_res_min                     = Sim_Struct.High_res_min;
gauss_filter_HighRes             = Sim_Struct.gauss_filter_HighRes;
larss_filter_HighRes             = Sim_Struct.larss_filter_HighRes;


if ~Sim_Struct.FORCE_SERIAL
    parfor i = 1:num_iterations
        for j = 1 : num_averages
            
            %% Filter in high resolution, then subsample
            
            High_2_Low_factor          = min_interval(i) / High_res_min;
            % Filter the delayed AIF with the gaussian kernel
            Temp_Ct_gauss_kernel       = filter(gauss_filter_HighRes(:,i)*High_res_min,1,Sim_AIF_HighRes_delayed_no_noise(:,i,j));
            % Filter the delayed AIF with Larsson's kernel
            Temp_Ct_larss_kernel       = filter(larss_filter_HighRes(:,i)*High_res_min,1,Sim_AIF_HighRes_delayed_no_noise(:,i,j));
            Sim_Ct_gauss_kernel(:,i,j) = downsample(Temp_Ct_gauss_kernel,High_2_Low_factor); %[mM]
            Sim_Ct_larss_kernel(:,i,j) = downsample(Temp_Ct_larss_kernel,High_2_Low_factor); %[mM]
            
            %         % Filter the delayed AIF with the gaussian kernel
            %         Sim_Ct_gauss_kernel(:,i,j) = filter(gauss_filter(:,i)*min_interval(i),1,Sim_AIF_delayed_no_noise(:,i,j));
            %         % Filter the delayed AIF with Larsson's kernel
            %         Sim_Ct_larss_kernel(:,i,j) = filter(larss_filter(:,i)*min_interval(i),1,Sim_AIF_delayed_no_noise(:,i,j));
            
            % Adding noise to simulated Ct(t) from both kernels
            noise_sigma_gauss         = mean(Sim_Ct_gauss_kernel(:,i,j)) ./ SNR_ratio(i);
            noise_to_add_gauss(:,i,j) = noise_sigma_gauss * randn(size(Sim_Ct_gauss_kernel(:,i,j)));
            noise_sigma_larss         = mean(Sim_Ct_larss_kernel(:,i,j))/SNR_ratio(i);
            noise_to_add_larss(:,i,j) = noise_sigma_larss * randn(size(Sim_Ct_larss_kernel(:,i,j)));
            
        end
        
        if ( mod(i,100) == 0 )
            display(sprintf('-I- Finished filtering 100 AIFs...'));
        end
        
        
    end
else
    for i = 1:num_iterations
        for j = 1 : num_averages
            
            %% Filter in high resolution, then subsample
            
            High_2_Low_factor          = min_interval(i) / High_res_min;
            % Filter the delayed AIF with the gaussian kernel
            Temp_Ct_gauss_kernel       = filter(gauss_filter_HighRes(:,i)*High_res_min,1,Sim_AIF_HighRes_delayed_no_noise(:,i,j));
            % Filter the delayed AIF with Larsson's kernel
            Temp_Ct_larss_kernel       = filter(larss_filter_HighRes(:,i)*High_res_min,1,Sim_AIF_HighRes_delayed_no_noise(:,i,j));
            Sim_Ct_gauss_kernel(:,i,j) = downsample(Temp_Ct_gauss_kernel,High_2_Low_factor); %[mM]
            Sim_Ct_larss_kernel(:,i,j) = downsample(Temp_Ct_larss_kernel,High_2_Low_factor); %[mM]
            
            %         % Filter the delayed AIF with the gaussian kernel
            %         Sim_Ct_gauss_kernel(:,i,j) = filter(gauss_filter(:,i)*min_interval(i),1,Sim_AIF_delayed_no_noise(:,i,j));
            %         % Filter the delayed AIF with Larsson's kernel
            %         Sim_Ct_larss_kernel(:,i,j) = filter(larss_filter(:,i)*min_interval(i),1,Sim_AIF_delayed_no_noise(:,i,j));
            
            % Adding noise to simulated Ct(t) from both kernels
            noise_sigma_gauss         = mean(Sim_Ct_gauss_kernel(:,i,j)) ./ SNR_ratio(i);
            noise_to_add_gauss(:,i,j) = noise_sigma_gauss * randn(size(Sim_Ct_gauss_kernel(:,i,j)));
            noise_sigma_larss         = mean(Sim_Ct_larss_kernel(:,i,j))/SNR_ratio(i);
            noise_to_add_larss(:,i,j) = noise_sigma_larss * randn(size(Sim_Ct_larss_kernel(:,i,j)));
            
        end
        
        if ( mod(i,100) == 0 )
            display(sprintf('-I- Finished filtering %d AIFs...',i));
        end
        
    end
end


end