function Out=ToShortName(in)
notcell=0;
if(~iscell(in))
    in={in};
    notcell=1;
end
for i=1:length(in)
%     S=regexp(in{i},'[ _]','split');
    S=regexpsplit(in{i},'[ _]');
    Out{i}=[upper(S{1}(1)) lower(S{1}(min(numel(S{1}),2))) upper(S{end}(1)) lower(S{end}(2))];
end
if(notcell)
    Out=Out{1};
end