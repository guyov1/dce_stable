function setComputerParamM(ParamName, Value)
%% Save CompParams in CompParams.mat
ParamName=lower(ParamName);
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
    CompParams=struct();
else
    CompParams=load('CompParams.mat');
    CompParams=CompParams.CompParams;
end
%%
CompParams.(ParamName)=Value;
save(CompParams_File_Path,'CompParams');
