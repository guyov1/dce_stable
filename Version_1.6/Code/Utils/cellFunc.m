function Out=cellFunc(Func,Vars)
N=numel(Vars(:));
Out=cell(N,1);
for i=1:N
    Out{i}=feval(Func,Vars{i});
end
Out=reshape(Out,size(Vars));