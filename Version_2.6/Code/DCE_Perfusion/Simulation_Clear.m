
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
clear structs;
clearvars -except Verbosity

load('myBreakpoints.mat');
delete('myBreakpoints.mat');
dbstop(current_break_points);

% Add current directory and all it's subdirectories to Matlab's path
addpath(genpath('./'));

% Close all open PDF files (so I will be able to write the report PDF later
if ~strcmp(Verbosity,'None')
    display('-I- Closing all open PDF files.');
    
end
[status, result] = dos('taskkill /F /IM acroRd32.exe');

% Prepare log file (pdf output file)
if (~exist('./Run_Output','dir'))
    mkdir('./Run_Output');
end
LogFN       = './Run_Output/Log.mat';
% Delete old file
if(exist(LogFN,'file'))
    delete(LogFN);
end
SN          = '\\underline{Report}';
Log.idx_000 = {['\\title{' SN '}\r\n\\maketitle\r\n']};
save(LogFN,'Log');
