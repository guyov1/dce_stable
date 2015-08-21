function Sim_Struct = setETMParams(Sim_Struct)
%setETMParams Summary of this function goes here

% Read the latest run parameters from the excel file
Sim_Struct.readLatestData                = true;
% Excel file to write and read from
Sim_Struct.ETM_filename                  = 'exported_kep_simulation.csv';
% Index to plot
Sim_Struct.ETM_idx_to_plot               = 14;
% Extended time to plot
Sim_Struct.total_sim_time_min_to_plot    = 20;
Sim_Struct.num_time_stamps_to_plot       = round(Sim_Struct.total_sim_time_min_to_plot / Sim_Struct.min_interval);
Sim_Struct.time_vec_minutes_to_plot      = (0 : Sim_Struct.num_time_stamps_to_plot - 1).* Sim_Struct.min_interval;

% Use a manual time vector for ETM estimation
Sim_Struct.use_manual_time_vec           = false;
if Sim_Struct.use_manual_time_vec
    Sim_Struct.manual_time_vec_minutes       = [0.1 0.2 0.4 0.5 1];
    Sim_Struct.time_vec_minutes              = Sim_Struct.manual_time_vec_minutes; % Adjust the simulation time vector
end
% Parameter just for ETM estimation
Sim_Struct.Ktrans_ETM_single = 0.11;
Sim_Struct.Ktrans_ETM_low    = 0.56; % 0.05
Sim_Struct.Ktrans_ETM_max    = 0.8; % 0.5
% Sim_Struct.Ktrans_ETM_vec    = linspace(Sim_Struct.Ktrans_ETM_low, Sim_Struct.Ktrans_ETM_max, Sim_Struct.num_iterations); % When iterating
Sim_Struct.Vp_ETM_single     = 0.1; 
Sim_Struct.Vp_ETM_low        = 0.1; % 0.01
Sim_Struct.Vp_ETM_max        = 0.2; % 0.3
% Sim_Struct.Vp_ETM_vec        = linspace(Sim_Struct.Vp_ETM_low, Sim_Struct.Vp_ETM_max, Sim_Struct.num_iterations); % When iterating

Sim_Struct.kep_ETM_single     = 0.01; % 0.01
Sim_Struct.kep_ETM_low        = 0.15; % 0.01
Sim_Struct.kep_ETM_max        = 0.25; % 0.3
% Sim_Struct.kep_ETM_vec        = linspace(Sim_Struct.kep_ETM_low, Sim_Struct.kep_ETM_low, Sim_Struct.num_iterations); % When iterating

%Sim_Struct.Ve_ETM_single     = 0.3;          
%Sim_Struct.Ve_ETM_low        = 0.01; % 0.01
%Sim_Struct.Ve_ETM_max        = 0.3; % 0.3
%Sim_Struct.Ve_ETM_vec        = linspace(Sim_Struct.Ve_ETM_low, Sim_Struct.Ve_ETM_low, Sim_Struct.num_iterations); % When iterating
Sim_Struct.constraint  = 'Rcheck = ( (Sim_Struct.Ve_ETM + Sim_Struct.Vp_ETM) > 1 ) | ( Sim_Struct.Ve_ETM > 1 ) | ( Sim_Struct.Vp_ETM > 1 );';
Sim_Struct.constraint  = 'Rcheck = ( Sim_Struct.Vp_ETM > 1 );';

if Sim_Struct.ETM_Model
   display('-I- ETM Ranges:'); 
   display([num2str(Sim_Struct.Ktrans_ETM_low) ' < Ktrans <' num2str(Sim_Struct.Ktrans_ETM_max)]);
   display([num2str(Sim_Struct.kep_ETM_low)    ' < kep    <' num2str(Sim_Struct.kep_ETM_max)]);
   display([num2str(Sim_Struct.Vp_ETM_low)     ' < Vp     <' num2str(Sim_Struct.Vp_ETM_max)]);
   % Ve derives from the rest
   ve_min = Sim_Struct.Ktrans_ETM_low / Sim_Struct.kep_ETM_max;
   ve_max = Sim_Struct.Ktrans_ETM_max / Sim_Struct.kep_ETM_low;
   display([num2str(ve_min) ' < Ve     <' num2str(ve_max)]);
end


end

