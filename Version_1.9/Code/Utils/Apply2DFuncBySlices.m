function Out=Apply2DFuncBySlices(In,Func)
Out=In;
for i=1:size(In,3)
    Out(:,:,i)=Func(squeeze(In(:,:,i)));
end