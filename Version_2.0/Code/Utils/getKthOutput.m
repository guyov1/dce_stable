function Out=getKthOutput(K,Func,Args)
[All{1:K}]=Func(Args{:});
Out=All{K};