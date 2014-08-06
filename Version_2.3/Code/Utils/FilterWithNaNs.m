function Out=FilterWithNaNs(I,Func,Thresh)
B=isfinite(I);
BS=Func(double(B));
I(~B)=0;
Out=Func(I)./BS;
Out(BS<Thresh)=NaN;