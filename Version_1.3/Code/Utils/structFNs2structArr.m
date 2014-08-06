function Out=structFNs2structArr(in)
FNs=fieldnames(in);
for i=1:length(FNs)
    Out(i)=in.(FNs{i});
end
