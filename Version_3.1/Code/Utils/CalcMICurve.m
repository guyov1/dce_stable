function MI=CalcNiiCurve(FNs,FNB)
MI=zeros(1,numel(FNs));
for i=1:numel(FNs)
    MI(i)=CalcNiiMI(FNs{i},FNB);
end