function Out=nCr(N,K)
dbif(N<171)
    Out=factorial(N)/(factorial(K)*factorial(N-K));
    return;
end
Out=1;
for i=N:-1:N-K+1
    Out=Out*i;
end
Out=Out/factorial(K);