function csvcwrite(C,FN)
fid=fopen(FN,'w');
for i=1:size(C,1)
    OutLine='';
    for j=1:size(C,2)
        if(ischar(C{i,j}))
            Cur=C{i,j};
        end
        if(isnumeric(C{i,j}))
            Cur=num2str(C{i,j});
        end
        OutLine=[OutLine Cur ','];
    end
    OutLine=OutLine(1:end-1);
    fprintf(fid,'%s\n',OutLine);
end
fclose(fid);
