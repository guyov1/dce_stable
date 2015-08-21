function Out=getDCMFile(in)
DD=dir([in filesep 'MR0*.dcm']);
if(isempty(DD))
    DD=dir([in filesep '00*.dcm']);
end
if(isempty(DD))
    DD=dir([in filesep 'MR1*.dcm']);
end
if(isempty(DD))
    Out=[];
else
    Out=[in filesep DD(2).name];
end