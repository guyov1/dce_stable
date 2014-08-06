% function [Cost More]=gCost(Raw,App,Method,Params)
% Method - 'SumAbs','RMS','SumLp','Lp','Corr','MAD','Mean'
% Compensation - 'AddCompensate' 'MultiCompensate' 'Relative'
function [Cost More k]=gCost(Raw,App,Method,Params)
if(~iscell(Method))
    Method={Method};
end
k=-777;
More=-777;
if(ismember('AddCompensate',Method))
    if(ismember('MultiCompensate',Method))
        k=[App App*0+1]\Raw;
        More=App*k(1)+k(2);
        AbsDiff=abs(Raw-More);
    else
        AbsDiff=abs(Raw-App-mean(Raw)+mean(App));
        k=mean(Raw)-mean(App);
        More=App+k;
    end
else
    if(ismember('MultiCompensate',Method))
        k=App\Raw;
        More=App*k;
        AbsDiff=abs(Raw-More);
    else
        AbsDiff=abs(Raw-App);
    end
end
if(ismember('Relative',Method))
    AbsDiff=AbsDiff./Raw;
end
MethodX=Method{ismember(Method,{'SumAbs','RMS','SumLp','Lp','Corr','MAD','Mean'})};
switch(MethodX)
    case 'SumAbs'
        Cost=sum(AbsDiff);
    case 'Mean'
        Cost=mean(AbsDiff);
    case 'RMS'
        Cost=sqrt(mean(AbsDiff.^2));
    case 'SumLp'
        Cost=sum(AbsDiff.^Params);
    case 'Lp'
        Cost=(sum(AbsDiff.^Params)).^(1/Params);
    case 'Corr'
        Cost=getKthElement(corrcoef(Raw,App),2);
    case 'MAD'
        Cost=median(AbsDiff);
end