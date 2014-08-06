function Out=getComputerParamM(ParamName, varargin)
CompParams=load('CompParams.mat');
CompParams=CompParams.CompParams;
Out=[];
ParamName=lower(ParamName);
if(isfield(CompParams,ParamName))
    Out=CompParams.(ParamName);
    if(~isempty(varargin))
        if(~isempty(varargin{1}))
            Out=getKthElement(Out,varargin{1});
        end
    end
else
        disp('Error getComputerParamsM');
        error('Error getComputerParamsM');
end