function Out=getKthElement(In,k)
if(iscell(In))
    Out=In{k};
    return;
else
    Out=In(k);
    return;
end
error('Error - getKthElement');