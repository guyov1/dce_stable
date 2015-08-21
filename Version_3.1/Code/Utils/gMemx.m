a=whos;
[S Ord]=sort([a.bytes],'descend');
for i=1:10
    disp([a(Ord(i)).name '    ' num2str(a(Ord(i)).bytes)]);
end
disp(['Total: ' num2str(sum(S)./(10^6)) 'MB']);