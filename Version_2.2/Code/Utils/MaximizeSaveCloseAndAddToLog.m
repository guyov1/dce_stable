function MaximizeSaveCloseAndAddToLog(WorkingP,FigName,LogHandle,LogTxt)
set(gcf,'Position',figposition([0 0 100 100]));
saveas(gcf,[WorkingP FigName '.png']);
saveas(gcf,[WorkingP FigName '.fig']);
close(gcf);
if(nargin>2)
    AddToLog(WorkingP,LogHandle,LogTxt,[FigName '.png']);
end