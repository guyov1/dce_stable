function In=UntilDelimiter(In,Del)
if(iscell(In))
    for i=1:numel(In)
        In{i}=UntilDelimiter(In{i},Del);
    end
else
    In=In(1:find(In==Del,1)-1);
end