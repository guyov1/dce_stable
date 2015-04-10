function gfig(a)
if(nargin>0)
    figure(a);
else
    figure;
end
set(gcf,'Position',figposition([10 10 70 70]));