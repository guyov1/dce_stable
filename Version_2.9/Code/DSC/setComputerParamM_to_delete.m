function setComputerParamM(ParamName, Value)
ParamName=lower(ParamName);
CompParams=load('CompParams.mat');
CompParams=CompParams.CompParams;
CompParams.(ParamName)=Value;

%% Save CompParams in CompParams.mat

% Get PC name
PC_name = getenv('ComputerName');

BaseBaseP=[pwd filesep];

% Get Matlab user name
TmpName= license('inuse');
TmpName=TmpName(1).user;

CompParams_File_Path =[BaseBaseP TmpName filesep 'Database' filesep 'CompParams.mat'];

save(CompParams_File_Path,'CompParams');



