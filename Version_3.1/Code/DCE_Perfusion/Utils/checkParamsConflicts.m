function [ Sim_Struct ] = checkParamsConflicts( Sim_Struct, Verbosity)
%checkParamsConflicts Check for parameters conflicts

% Avoid memory overhead if there are too many iteratins
if (Sim_Struct.num_iterations > 40)
    Sim_Struct.FORCE_SERIAL                  = true;
    Sim_Struct.FORCE_MAIN_LOOP_SERIAL        = true;
end

% Ignoring gaussian calculation but wanting to iterate over it's parameters
if (Sim_Struct.Ignore_Gaussian_Calculation && (Sim_Struct.iterate_gaussian_sigma || Sim_Struct.iterate_gaussian_time_delay || Sim_Struct.iterate_gaussian_time_delay) )
    error('Set to ignore Gaussian calculation. Cant iterate over its parameters!');
end

% One-hot option only for de-convolution method
if (Sim_Struct.Use_Cyclic_Conv_4_ht_est + Sim_Struct.Use_Upsampling_Delay_Comp + Sim_Struct.Correct_estimation_due_to_delay > 1)
    error('More than 1 option for deconvolution is set!');
end

% Usampling + Cyclic can only be set in case Cyclic De-convolve is set
if (Sim_Struct.Use_Upsampling_and_Cyclic && ~Sim_Struct.Use_Cyclic_Conv_4_ht_est)
    error('Cyclic De-convolve is unset! Cant up-sample!');
end

% Use AIF delay corrcetion only in case we allowed it
if (Sim_Struct.LQ_Model_AIF_Delay_Correct && ~Sim_Struct.Correct_estimation_due_to_delay)
    error('Cant correct for AIF delay, because correction is not enabled!');
end

% Use AIF delay corrcetion only in case we allowed it
if (Sim_Struct.Simple_AIF_Delay_Correct && ~Sim_Struct.Correct_estimation_due_to_delay)
    error('Cant correct for AIF delay, because correction is not enabled!');
end

% Use AIF delay corrcetion only in case we allowed it
if (Sim_Struct.Simple_AIF_Delay_Correct + Sim_Struct.LQ_Model_AIF_Delay_Correct > 1)
    error('More than one AIF delay correction method chosen!');
end

if (Sim_Struct.ETM_Model && Sim_Struct.Adjusted_Larsson_Model) || (~Sim_Struct.ETM_Model && ~Sim_Struct.Adjusted_Larsson_Model)
    error('Choose either ETM or Adjusted Larsson Model!');
end

if strcmp(Verbosity,'Full')
    display('-I- Finished Setting simulation parameters...');
end

end

