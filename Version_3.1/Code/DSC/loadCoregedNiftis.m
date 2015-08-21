function [Out Headers]=loadCoregedNiftis(FNs,Xs,Ys,Zs)
[Out(:,:,:,1) Headers{1}]=loadniidata(FNs{1});
S=size(Out);
Partial=nargin>1;
n=length(FNs);
if(n>1)
    Headers{n}=[];
end
if(Partial)
    if(~exist('Xs','var'))
        Xs=1:S(1);
    else
        if(isempty(Xs))
            Xs=1:S(1);
        end
    end
    if(~exist('Ys','var'))
        Ys=1:S(2);
    else
        if(isempty(Ys))
            Ys=1:S(2);
        end
    end
    if(~exist('Zs','var'))
        Zs=1:S(3);
    else
        if(isempty(Zs))
            Zs=1:S(3);
        end
    end
    Out=Out(Xs,Ys,Zs);
    Out(end,end,end,n)=0;
    for i=2:n
%         disp(i);
        [Tmp Headers{i}]=loadniidata(FNs{i});
        Out(:,:,:,i)=Tmp(Xs,Ys,Zs);
    end
else
    Out(end,end,end,n)=0;
    for i=2:n
%         disp(i);
        [Out(:,:,:,i) Headers{i}]=loadniidata(FNs{i});
    end
end