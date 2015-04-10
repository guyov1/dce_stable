function In=UntilDelimiter(In,Del,K)
if(nargin<3)
    K=1;
end
if(iscell(In))
    for i=1:numel(In)
        In{i}=UntilDelimiter(In{i},Del,K);
    end
else
    if(K>0)
        In=In(1:find(In==Del,K)-1);
    else
        In=In(1:find(In==Del,abs(K),'last')-1);
    end
end