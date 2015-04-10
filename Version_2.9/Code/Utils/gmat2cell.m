function C=gmat2cell(M,Method)
if(~exist('Method','var'))
    Method=0;
end
switch Method
    case 0
        C=mat2cell(M, ones(1,size(M,1)), ones(1,size(M,2)));
    case 1
        C=mat2cell(M, ones(1,size(M,1)), size(M,2));
end