function Out=num2strC(in,precision)
Out=cell(size(in));
for i=1:numel(in)
    if(nargin>1)
        Out{i}=num2str(in(i),precision);
    else
        Out{i}=num2str(in(i));
    end
end