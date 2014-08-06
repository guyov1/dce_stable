function Out=CorrectStrForLog(In)

In=strrep(In,'_','-');
Out=strrep(In,'\','/');