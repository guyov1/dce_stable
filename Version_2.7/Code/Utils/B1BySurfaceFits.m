SurfFit=cell(1,size(RefVolX,3));
for i=1:size(RefVolX,3)
    disp(i);
    F=find(I{3}==i);
    if(isempty(F))
        continue;
    end
    SurfFit{i}=fit([I{1}(F) I{2}(F)],SolVec(F), 'poly11' );
end
%%
for i=1:size(RefVolX,3)
    disp(i);
    F=find(IF{3}==i);
    if(isempty(SurfFit{i}))
        continue;
    end
    Tmp=squeeze(NewField(:,:,i));
    Tmp(:)=feval(SurfFit{i},[IF{1}(F) IF{2}(F)]);
    NewField(:,:,i)=Tmp;
end
CT1=(UncleanedT1./NewField)*RefTrgVal;