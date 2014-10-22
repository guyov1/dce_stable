function Out=regexpsplit(Str,Exp)
Len=numel(Str);
[Start End]=regexp(Str,Exp);
NStart=[1 (End+1)];
NEnd=[(Start-1) Len];
NOccur=numel(NStart);
Out=cell(1,NOccur);
for i=1:NOccur
    Out{i}=Str(NStart(i):NEnd(i));
end
Out=Out(~gIsEmpty(Out));