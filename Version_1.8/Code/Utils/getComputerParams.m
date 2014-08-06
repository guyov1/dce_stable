function Out=getComputerParams(ParamName, varargin)

CompParams=load('CompParams.mat');
CompParams=CompParams.CompParams;
ParamName=lower(ParamName);

if(isfield(CompParams,ParamName))
    Out=CompParams.(ParamName);
    if(~isempty(varargin))
        if(~isempty(varargin{1}))
            Out=getKthElement(Out,varargin{1});
        end
    end
else
    Out=[];
    disp('Error getComputerParams');
    error('Error getComputerParams');
end

