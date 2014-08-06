function Out=repMulti(A,B)

if(all(size(A)==size(B)))
    Out=A.*B;
    return;
end

% We assume that either the lines or the columns are equal
if(size(A,1)==size(B,1))
    Out=A.*repmat(B,1,size(A,2));
else
    Out=A.*repmat(B,size(A,1),1);
end