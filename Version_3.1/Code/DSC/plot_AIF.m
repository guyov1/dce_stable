function [] = plot_AIF(AIF,CTC,handles)

axes(handles.plot_AIFs);

% choose the max of y axis, for a comfortable view (there are 2
% different options. One is the global maximum, the other is the max of
% ~90% from the curves.
local_max=max(CTC);
if local_max>handles.c_t_majority_max
    y_axis_max=handles.c_t_global_max;
else
    y_axis_max=handles.c_t_majority_max;
end
% plot(AIF,'r','Linewidth',1.5);
% axis([handles.first_bl_sample handles.last_sample handles.c_t_global_min y_axis_max]);

plot(AIF,'r','LineWidth',1.5);
    hold on;
    plot(CTC);
    hold off;
    grid on
    axis([handles.first_bl_sample handles.last_sample handles.c_t_global_min y_axis_max]);
    set(gca,'xtick',[handles.first_bl_sample:handles.last_bl_sample-handles.first_bl_sample+1:handles.last_bl_sample]);
