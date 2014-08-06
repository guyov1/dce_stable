function setComputerParamM(ParamName, Value)
CompParams=load('CompParams.mat');
CompParams=CompParams.CompParams;
CompParams.(ParamName)=Value;
save('CompParams.mat','CompParams');