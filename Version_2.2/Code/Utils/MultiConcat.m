function Out=MultiConcat(A,B)
if(iscell(B))
    for i=1:length(B)
        Out{i}=[A B{i}];
    end
else
    for i=1:length(A)
        Out{i}=[A{i} B];
    end
end