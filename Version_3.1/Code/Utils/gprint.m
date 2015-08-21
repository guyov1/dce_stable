% function gprint(Handle,FN,R)
function gprint(Handle,FN,R)
% set(gcf,'Position',figposition([0 0 100 100]));
if(~exist('R','var'))
    R=[];
end
if(isempty(R))
    R=150;
    print(FN,['-f' num2str(Handle)],'-dpng');    
else
    print(FN,['-f' num2str(Handle)],'-djpeg',['-r' num2str(R)]);
end
