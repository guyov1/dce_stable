function setComputerParamM(ParamName, Value)
ParamName=lower(ParamName);
CompParams=load('CompParams.mat');
CompParams=CompParams.CompParams;
CompParams.(ParamName)=Value;

%% Save CompParams in CompParams.mat
CompParams_File_Path=which('CompParams.mat');
if(isempty(CompParams_File_Path))
    % Get PC name
    PC_name = getenv('ComputerName');
    
    BaseBaseP=[pwd filesep];
    
    % Get Matlab user name
    %TmpName= license('inuse');
    %TmpName=TmpName(1).user;
    TmpName = char(java.lang.System.getProperty('user.name'));
    
    CompParams_File_Path =[BaseBaseP TmpName filesep 'Database' filesep 'CompParams.mat'];
end
save(CompParams_File_Path,'CompParams');


