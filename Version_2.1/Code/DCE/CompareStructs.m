function Out=CompareStructs(A,B)
% A=FullInfos(1);
% B=FullInfos(3);
FldNmsA=fieldnames(A);
FldNmsB=fieldnames(B);
FldDiff=union(setdiff(FldNmsA,FldNmsB),setdiff(FldNmsB,FldNmsA));
IFlds=intersect(FldNmsA,FldNmsB);
%%
Out=cell(0,1);
for i=1:numel(IFlds)
    ValA=A.(IFlds{i});
    ValB=B.(IFlds{i});
    if(ischar(ValA))
        if(~strcmp(ValA,ValB))
            disp([IFlds{i} ' ' ValA ' ' ValB]);
            Out{end+1}=IFlds{i};
        end
    else
        if(isstruct(ValA))
            CompareStructs(ValA,ValB);
        else
            if(all(size(ValA)==size(ValB)))
                if(any(ValA~=ValB))
                    disp([IFlds{i} ' se' num2str(size(A)) ' ' num2str(ValA(1)) ' : ' num2str(ValB(1))]);
                    Out{end+1}=IFlds{i};
                end
            else
                disp([IFlds{i} ' sX' num2str(size(A)) ' ' 's' num2str(size(A)) ' ' num2str(ValA(1)) ' : ' num2str(ValB(1))]);
                Out{end+1}=IFlds{i};
            end
        end
    end
end