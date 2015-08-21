function [ ] = plotSingleCTC( x, y, z, CTC_4D, conv_result_IRF_4D)
%plotSingleCTC Plot a single CTC of the brain

figure;
hold on;
plot(squeeze(CTC_4D(x,y,z,:)),'b');
plot(squeeze(conv_result_IRF_4D(x,y,z,:)),'g');
hold off;

figure;
hold on;
h1 = plot(time_vec_minutes,squeeze(CTC_4D(x,y,z,:)),'LineWidth',6,'Color','k');
%h2 = plot(time_vec_minutes,AIF_part,'LineWidth',1,'LineStyle','+','Color','r');
%h3 = plot(time_vec_minutes,Kep_Filter_Part,'LineWidth',1,'LineStyle','o','Color','g');
%h4 = plot(time_vec_minutes,Sum_Result,'LineWidth',2,'Color','b');
%h5 = plot(time_vec_minutes,Sum_Result_NonLinear,'LineWidth',2,'Color','c');
hold off;
title('Concentration Time Curve','fontsize',15,'FontWeight','bold');
%legend([h1 h2 h3 h4],'CTC','Tofts AIF','Tofts permeability','Tofts fit')
xlabel('Time [Min]','fontsize',15,'FontWeight','bold');
ylabel('C_t(t) [mM]','fontsize',15,'FontWeight','bold');
set(gca,'fontsize',15,'FontWeight','bold');


end

