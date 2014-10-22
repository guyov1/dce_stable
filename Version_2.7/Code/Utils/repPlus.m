function Out=repPlus(A,B)
if(size(A,1)==size(B,1))
    Out=A+repmat(B,1,size(A,2));
else
    Out=A+repmat(B,size(A,1),1);
end