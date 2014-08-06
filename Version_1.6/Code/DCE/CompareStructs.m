A=FullInfos(1);
B=FullInfos(3);
FldNmsA=fieldnames(A);
FldNmsB=fieldnames(A);
FldDiff=union(setdiff(FldNmsA,FldNmsB),setdiff(FldNmsB,FldNmsA));
IFlds=intersect(FldNmsA,FldNmsB);
%%
for i=1:numel(IFlds)
    ValA=A.(IFlds{i});
    ValB=B.(IFlds{i});
    if(ischar(ValA))
        if(~strcmp(ValA,ValB))
            disp([IFlds{i} ' ' ValA ' ' ValB]);
        end
    else
        if(isstruct(ValA))
        else
            if(any(ValA~=ValB))
                if(all(size(A)==size(B)))
                    disp([IFlds{i} ' se' num2str(size(A)) ' ' num2str(ValA(1)) ' : ' num2str(ValB(1))]);
                else
                    disp([IFlds{i} ' sX' num2str(size(A)) ' ' 's' num2str(size(A)) ' ' num2str(ValA(1)) ' : ' num2str(ValB(1))]);
                end
            end
        end
    end
end