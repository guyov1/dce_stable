function Out=flattenCell(in,Dim)
% in=CrgFN;
S=size(in);
if(~exist('Dim','var'))
    Dim=find(S~=1);
end

BC=cellFuncX('iscell',in);
CN=cellNumel(in);
CCN=CN;
CCN(BC==0)=1;
Total=sum(CCN);
nS=[S(1:Dim-1) Total S(Dim+1:end)];
Out=cell(nS);
k=1;
for i=1:numel(CCN)
    if(BC(i))
        Out(k:(k+CCN(i)-1))=in{i};
    else
        Out(k)=in(i);
    end
    k=k+CCN(i);
end