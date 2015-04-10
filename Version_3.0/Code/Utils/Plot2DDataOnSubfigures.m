function Plot2DDataOnSubfigures(Handle,Data)

figure(Handle);clf;
for i=1:size(Data,1)
    gsubplot(size(Data,1),i);
    plot(Data(i,:));
end

