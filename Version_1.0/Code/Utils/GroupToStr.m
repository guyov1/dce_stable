function Out=GroupToStr(in)
Out='[';
for i=1:numel(in)-1
    Out=[Out in{i} '_'];
end
Out=[Out in{end} ']'];