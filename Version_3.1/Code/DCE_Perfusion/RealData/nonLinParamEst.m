function [ ...
    Flow_with_Delay_vec, Flow_no_Delay_vec, Delay_sec_by_Max_Val_with_Delay, Delay_sec_by_Max_Val_no_Delay, ...
    fitted_larsson_with_Delay, fitted_larsson_no_Delay, fitted_larsson_with_Delay_High_F, fitted_larsson_no_Delay_High_F, fitted_larsson_with_Delay_no_Ve, fitted_larsson_no_Delay_no_Ve, ...
    fitted_larsson_with_Delay_no_E, fitted_larsson_no_Delay_no_E, fitted_larsson_with_Delay_no_E_High_F, fitted_larsson_no_Delay_no_E_High_F, ...
    fitted_gaussian, fitted_double_gaussian, gaussian_param_vec, double_gaussian_param_vec, ...
    Ktrans_with_Delay_vec, Ktrans_no_Delay_vec, ...
    E_with_Delay_vec, E_no_Delay_vec, E_with_Delay_no_Ve_vec, E_no_Delay_no_Ve_vec, Ktrans_with_Delay_High_F_vec, Ktrans_no_Delay_High_F_vec, ...
    Vb_with_Delay_vec, Vb_no_Delay_vec, Vb_with_Delay_High_F_vec, Vb_no_Delay_High_F_vec, Vb_with_Delay_no_Ve_vec, Vb_no_Delay_no_Ve_vec,...
    Vb_with_Delay_no_E_vec, Vb_no_Delay_no_E_vec, Vb_with_Delay_no_E_High_F_vec, Vb_no_Delay_no_E_High_F_vec, ...
    Ve_with_Delay_vec, Ve_no_Delay_vec, Ve_with_Delay_High_F_vec, Ve_no_Delay_High_F_vec, MTT_with_Delay_vec, MTT_no_Delay_vec, ...
    Ktrans_with_Delay_Patlak_vec, Ktrans_no_Delay_Patlak_vec, Vb_with_Delay_Patlak_vec, Vb_no_Delay_Patlak_vec, ...
    MTT_with_Delay_Patlak_vec, MTT_no_Delay_Patlak_vec ] ...
    = nonLinParamEst(Sim_Struct,Est_ht_with_Delay,Est_ht_no_Delay,Ct,AIF_delay_corrected,AIF_no_Delay,idx_fig, Parallel_Real_Data_Est)

display('--------------------------------------------------------');
display('-I- Starting non linear parameters estimation...');
display('--------------------------------------------------------');


num_voxels                            = Sim_Struct.num_voxels;
num_time_stamps                       = Sim_Struct.num_time_stamps;

%% Initiate return vectors
% Gaussian parameters
gaussian_param_vec                    = zeros(3,num_voxels);
double_gaussian_param_vec             = zeros(6,num_voxels);

fitted_larsson_with_Delay             = zeros(num_voxels,num_time_stamps);
fitted_larsson_no_Delay               = zeros(num_voxels,num_time_stamps);
fitted_larsson_with_Delay_High_F      = zeros(num_voxels,num_time_stamps);
fitted_larsson_no_Delay_High_F        = zeros(num_voxels,num_time_stamps);
fitted_larsson_with_Delay_no_Ve       = zeros(num_voxels,num_time_stamps);
fitted_larsson_no_Delay_no_Ve         = zeros(num_voxels,num_time_stamps);
fitted_larsson_with_Delay_no_E        = zeros(num_voxels,num_time_stamps);
fitted_larsson_no_Delay_no_E          = zeros(num_voxels,num_time_stamps);
fitted_larsson_with_Delay_no_E_High_F = zeros(num_voxels,num_time_stamps);
fitted_larsson_no_Delay_no_E_High_F   = zeros(num_voxels,num_time_stamps);

fitted_gaussian                       = zeros(num_voxels,num_time_stamps);
fitted_double_gaussian                = zeros(num_voxels,num_time_stamps);

% Larsson parameters
Flow_with_Delay_vec                   = zeros(1,num_voxels);
Flow_no_Delay_vec                     = zeros(1,num_voxels);

Delay_sec_by_Max_Val_with_Delay       = zeros(1,num_voxels);
Delay_sec_by_Max_Val_no_Delay         = zeros(1,num_voxels);

Ktrans_with_Delay_vec                 = zeros(1,num_voxels);
Ktrans_no_Delay_vec                   = zeros(1,num_voxels);
Ktrans_with_Delay_Patlak_vec          = zeros(1,num_voxels);
Ktrans_no_Delay_Patlak_vec            = zeros(1,num_voxels);
Ktrans_with_Delay_High_F_vec          = zeros(1,num_voxels);
Ktrans_no_Delay_High_F_vec            = zeros(1,num_voxels);

E_with_Delay_vec                      = zeros(1,num_voxels);
E_no_Delay_vec                        = zeros(1,num_voxels);
E_with_Delay_no_Ve_vec                = zeros(1,num_voxels);
E_no_Delay_no_Ve_vec                  = zeros(1,num_voxels);

Vb_with_Delay_vec                     = zeros(1,num_voxels);
Vb_no_Delay_vec                       = zeros(1,num_voxels);

Vb_with_Delay_High_F_vec              = zeros(1,num_voxels);
Vb_no_Delay_High_F_vec                = zeros(1,num_voxels);
Vb_with_Delay_no_Ve_vec               = zeros(1,num_voxels);
Vb_no_Delay_no_Ve_vec                 = zeros(1,num_voxels);

Vb_with_Delay_no_E_vec                = zeros(1,num_voxels);
Vb_no_Delay_no_E_vec                  = zeros(1,num_voxels);
Vb_with_Delay_no_E_High_F_vec         = zeros(1,num_voxels);
Vb_no_Delay_no_E_High_F_vec           = zeros(1,num_voxels);
Vb_with_Delay_Patlak_vec              = zeros(1,num_voxels);
Vb_no_Delay_Patlak_vec                = zeros(1,num_voxels);

Ve_with_Delay_vec                     = zeros(1,num_voxels);
Ve_no_Delay_vec                       = zeros(1,num_voxels);
Ve_with_Delay_High_F_vec              = zeros(1,num_voxels);
Ve_no_Delay_High_F_vec                = zeros(1,num_voxels);

MTT_with_Delay_vec                    = zeros(1,num_voxels);
MTT_no_Delay_vec                      = zeros(1,num_voxels);
MTT_with_Delay_Patlak_vec             = zeros(1,num_voxels);
MTT_no_Delay_Patlak_vec               = zeros(1,num_voxels);

%% Estimation

if Parallel_Real_Data_Est
    parfor j=1:num_voxels
        [gaussian_param_vec(:,j), fitted_gaussian(j,:), double_gaussian_param_vec(:,j), fitted_double_gaussian(j,:),...
            Flow_with_Delay_vec(j), Flow_no_Delay_vec(j), Vb_with_Delay_vec(j), Vb_no_Delay_vec(j), E_with_Delay_vec(j), E_no_Delay_vec(j), ...
            Ktrans_with_Delay_vec(j), Ktrans_no_Delay_vec(j), Ve_with_Delay_vec(j), Ve_no_Delay_vec(j), MTT_with_Delay_vec(j), MTT_no_Delay_vec(j),...
            fitted_larsson_with_Delay(j,:), fitted_larsson_no_Delay(j,:),...
            Vb_with_Delay_High_F_vec(j), Ktrans_with_Delay_High_F_vec(j), Ve_with_Delay_High_F_vec(j), fitted_larsson_with_Delay_High_F(j,:),...
            Vb_no_Delay_High_F_vec(j), Ktrans_no_Delay_High_F_vec(j), Ve_no_Delay_High_F_vec(j), fitted_larsson_no_Delay_High_F(j,:),...
            Vb_with_Delay_no_Ve_vec(j), E_with_Delay_no_Ve_vec(j), fitted_larsson_with_Delay_no_Ve(j,:),...
            Vb_no_Delay_no_Ve_vec(j), E_no_Delay_no_Ve_vec(j), fitted_larsson_no_Delay_no_Ve(j,:),...
            Vb_with_Delay_no_E_vec(j), fitted_larsson_with_Delay_no_E(j,:),...
            Vb_no_Delay_no_E_vec(j), fitted_larsson_no_Delay_no_E(j,:),...
            Vb_with_Delay_no_E_High_F_vec(j), fitted_larsson_with_Delay_no_E_High_F(j,:),...
            Vb_no_Delay_no_E_High_F_vec(j), fitted_larsson_no_Delay_no_E_High_F(j,:),...
            Ktrans_with_Delay_Patlak_vec(j), Ktrans_no_Delay_Patlak_vec(j), Vb_with_Delay_Patlak_vec(j), Vb_no_Delay_Patlak_vec(j),...
            MTT_with_Delay_Patlak_vec(j), MTT_no_Delay_Patlak_vec(j), Delay_sec_by_Max_Val_with_Delay(j),Delay_sec_by_Max_Val_no_Delay(j)]...
            = Params_Est_Real_Data(Sim_Struct, Est_ht_with_Delay(:,j), Est_ht_no_Delay(:,j), Ct(j,:), AIF_delay_corrected(:,j), AIF_no_Delay, idx_fig);
        % Report after each 100 voxels
        if ( mod(j,1000) == 0 )
            display(sprintf('Finished lsqcurvefit for %d voxels...',j));
            remaining_voxels = num_voxels - j;
            fprintf('Number of remaining voxels: %d .\n',remaining_voxels);
        end
    end
else
    for j=1:num_voxels
        [gaussian_param_vec(:,j), fitted_gaussian(j,:), double_gaussian_param_vec(:,j), fitted_double_gaussian(j,:),...
            Flow_with_Delay_vec(j), Flow_no_Delay_vec(j), Vb_with_Delay_vec(j), Vb_no_Delay_vec(j), E_with_Delay_vec(j), E_no_Delay_vec(j), ...
            Ktrans_with_Delay_vec(j), Ktrans_no_Delay_vec(j), Ve_with_Delay_vec(j), Ve_no_Delay_vec(j), MTT_with_Delay_vec(j), MTT_no_Delay_vec(j),...
            fitted_larsson_with_Delay(j,:), fitted_larsson_no_Delay(j,:),...
            Vb_with_Delay_High_F_vec(j), Ktrans_with_Delay_High_F_vec(j), Ve_with_Delay_High_F_vec(j), fitted_larsson_with_Delay_High_F(j,:),...
            Vb_no_Delay_High_F_vec(j), Ktrans_no_Delay_High_F_vec(j), Ve_no_Delay_High_F_vec(j), fitted_larsson_no_Delay_High_F(j,:),...
            Vb_with_Delay_no_Ve_vec(j), E_with_Delay_no_Ve_vec(j), fitted_larsson_with_Delay_no_Ve(j,:),...
            Vb_no_Delay_no_Ve_vec(j), E_no_Delay_no_Ve_vec(j), fitted_larsson_no_Delay_no_Ve(j,:),...
            Vb_with_Delay_no_E_vec(j), fitted_larsson_with_Delay_no_E(j,:),...
            Vb_no_Delay_no_E_vec(j), fitted_larsson_no_Delay_no_E(j,:),...
            Vb_with_Delay_no_E_High_F_vec(j), fitted_larsson_with_Delay_no_E_High_F(j,:),...
            Vb_no_Delay_no_E_High_F_vec(j), fitted_larsson_no_Delay_no_E_High_F(j,:),...
            Ktrans_with_Delay_Patlak_vec(j), Ktrans_no_Delay_Patlak_vec(j), Vb_with_Delay_Patlak_vec(j), Vb_no_Delay_Patlak_vec(j),...
            MTT_with_Delay_Patlak_vec(j), MTT_no_Delay_Patlak_vec(j), Delay_sec_by_Max_Val_with_Delay(j),Delay_sec_by_Max_Val_no_Delay(j)]...
            = Params_Est_Real_Data(Sim_Struct, Est_ht_with_Delay(:,j), Est_ht_no_Delay(:,j), Ct(j,:), AIF_delay_corrected(:,j), AIF_no_Delay, idx_fig);
        % Report after each 100 voxels
        if ( mod(j,1000) == 0 )
            display(sprintf('Finished lsqcurvefit for %d voxels...',j));
            remaining_voxels = num_voxels - j;
            fprintf('Number of remaining voxels: %d .\n',remaining_voxels);
        end
    end
end


end


