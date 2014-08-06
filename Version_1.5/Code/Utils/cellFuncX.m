function Out=cellFuncX(Func,Vars)
Out=cellFunc(Func,Vars);
ZL=cellFunc('islogical',Out);
ZN=cellFunc('isnumeric',Out);
ZX=cellFunc('numel',Out);
if((all([ZL{:}]) || all([ZN{:}])) && all([ZX{:}]==1))
    Out=[Out{:}];
    Out=reshape(Out,size(Vars));
end