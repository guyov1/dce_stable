function Out=getDCMFile(in)
DD=dir([in filesep 'MR0*.dcm']);
if(isempty(DD))
    Out=[];
else
    Out=[in filesep DD(1).name];
end