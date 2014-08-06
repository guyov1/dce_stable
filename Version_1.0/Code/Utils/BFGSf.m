function [X,FVAL] = BFGSf(FUN,X0,LB,UB,optionsB,varargin)

ConstX=@(x,l,h) min(h,max(l,x));
BCostFunc=@(x) FUN(ConstX(x,LB,UB));
if(isempty(optionsB))
    optionsB = struct('GradObj','off','Display','iter');
end
[X,FVAL] = fminlbfgs(BCostFunc,X0,optionsB);
X=ConstX(X,LB,UB);