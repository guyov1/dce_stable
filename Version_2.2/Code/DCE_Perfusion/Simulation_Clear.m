
% Close all open figures and clear everything except breakpoints
if ~strcmp(Verbosity,'None')
    display('-I- Closing figures and cleaning...');
end

close all;
current_break_points = dbstatus;
save('myBreakpoints.mat', 'current_break_points');
% Clearing
%clear all;

clear functions;
clear global;
clear import;
clear mex;
clearvars -except Verbosity

load('myBreakpoints.mat');
delete('myBreakpoints.mat');
dbstop(current_break_points);

% Add current directory and all it's subdirectories to Matlab's path
addpath(genpath('./'));

% Enable parallel processing
try
    num_processes = getenv('NUMBER_OF_PROCESSORS');
    % Maximum allowed processes are 12 in matlab 2012
    if (num_processes > 12)
        num_processes = 12;
    end
    
    myCluster = parcluster('local');
    myCluster.NumWorkers = num_processes; % 'Modified' property now TRUE
    saveProfile(myCluster);   % 'local' profile now updated,
    
    if ~strcmp(Verbosity,'None')
        fprintf('\n');
        display('-I- Initiating matlab pool for parallel processing...');
        fprintf('\n');
    end
    
    matlabpool;
    
    if ~strcmp(Verbosity,'None')
        fprintf('\n');
        display('-I- Finished matlab pool initiation.');
        fprintf('\n');
        
    end
    
catch error
    
    if ~strcmp(Verbosity,'None')
        fprintf('\n');
        fprintf('\n');
        display('-I- Matlab pool already running!');
        fprintf('\n');
    end
    
end

% Close all open PDF files (so I will be able to write the report PDF later
if ~strcmp(Verbosity,'None')
    display('-I- Closing all open PDF files.');
    
end
[status, result] = dos('taskkill /F /IM acroRd32.exe');

% Prepare log file (pdf output file)
if (~exist('./Run_Output','dir'))
    mkdir('./Run_Output');
end
LogFN       = 'Run_Output/Log.mat';
% Delete old file
if(exist(LogFN,'file'))
    delete(LogFN);
end
SN          = '\\underline{Report}';
Log.idx_000 = {['\\title{' SN '}\r\n\\maketitle\r\n']};
save(LogFN,'Log');
