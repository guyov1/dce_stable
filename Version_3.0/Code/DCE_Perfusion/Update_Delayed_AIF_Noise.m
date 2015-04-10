function [ noise_to_add_delayed_AIF, noise_to_add_delayed_AIF_HighRes] = Update_Delayed_AIF_Noise( noise_to_add_AIF, noise_to_add_AIF_HighRes, delay_index, FORCE_SERIAL, num_iterations, num_averages )


%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

noise_to_add_delayed_AIF         = zeros(size(noise_to_add_AIF));
noise_to_add_delayed_AIF_HighRes = zeros(size(noise_to_add_AIF_HighRes));

if ~FORCE_SERIAL
    
    parfor iter_num = 1 : num_iterations
        for avg_num = 1 : num_averages
            
            % Right shift the AIF
            if (delay_index > 0)
                noise_to_add_delayed_AIF(:,iter_num,avg_num)         = Update_Pos_Delay(delay_index(iter_num), noise_to_add_AIF(:,iter_num,avg_num) );
                noise_to_add_delayed_AIF_HighRes(:,iter_num,avg_num) = Update_Pos_Delay(delay_index(iter_num), noise_to_add_AIF_HighRes(:,iter_num,avg_num) );
            elseif (delay_index < 0) % Left shift
                noise_to_add_delayed_AIF(:,iter_num,avg_num)         = Update_Neg_Delay(delay_index(iter_num), noise_to_add_AIF(:,iter_num,avg_num) );
                noise_to_add_delayed_AIF_HighRes(:,iter_num,avg_num) = Update_Neg_Delay(delay_index(iter_num), noise_to_add_AIF_HighRes(:,iter_num,avg_num) );
            else % no shift
                
            end
            
        end
    end 
    
else
    
    for iter_num = 1:num_iterations
        for avg_num = 1 : num_averages
            
            % Right shift the AIF
            if (delay_index > 0)
                noise_to_add_delayed_AIF(:,iter_num,avg_num)         = Update_Pos_Delay(delay_index(iter_num), noise_to_add_AIF(:,iter_num,avg_num) );
                noise_to_add_delayed_AIF_HighRes(:,iter_num,avg_num) = Update_Pos_Delay(delay_index(iter_num), noise_to_add_AIF_HighRes(:,iter_num,avg_num) );
            elseif (delay_index < 0) % Left shift
                noise_to_add_delayed_AIF(:,iter_num,avg_num)         = Update_Neg_Delay(delay_index(iter_num), noise_to_add_AIF(:,iter_num,avg_num) );
                noise_to_add_delayed_AIF_HighRes(:,iter_num,avg_num) = Update_Neg_Delay(delay_index(iter_num), noise_to_add_AIF_HighRes(:,iter_num,avg_num) );
            else % no shift
                
            end
            
        end
    end 
    
    
end

end

