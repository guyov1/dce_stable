BasePSub='\\fmri-t9\users\Moran\DCE\DataBase\SmVl_20120930\';
BasePFull='\\fmri-t9\users\Moran\DCE\DataBase\GB-HTR\SMOLIAR_VLADIMIR\SmVl_20120930\';
FA=loadniidata([BasePFull 'AutoArtBAT\F1v0.nii']);
FB=loadniidata([BasePSub 'AutoArtBAT\F1v0.nii']);
RMS=loadniidata([BasePFull 'AutoArtBAT\RSS0.nii']);
FThresh=0.990;
FMsk=FA>FThresh & FB>FThresh & RMS>0;
figure;imagesc(mritransform(FMsk(:,:,2)));
%%
BA=loadniidata([BasePFull 'AutoArtBAT\BATfinalSec.nii']);
BB=loadniidata([BasePSub 'AutoArtBAT\BATfinalSec.nii']);

Edg=unique(BA(FMsk));
HH=hist3([BA(FMsk) BB(FMsk)],{Edg,Edg});
figure(501);clf;
set(gca,'FontSize',20);
axis equal
imagesc(HH);
[R P]=corrcoef([BA(FMsk) BB(FMsk)]);
% title(getKthElement(corrcoef([BA(FMsk) BB(FMsk)]),2));
title(['HTR vs. STR 2D histogram, correlation=' num2str(R(2),'%2.2f'), ', p<1e-10']);
colormap gray
Tick=floor(linspace(1,ceil(numel(Edg)/2),3));
Tick=[Tick numel(Edg)-Tick(2)+1 numel(Edg)];
set(gca,'XTick',Tick);
set(gca,'YTick',Tick);
set(gca,'XTickLabel',num2strC(Edg(Tick),'%2.0f'));
set(gca,'YTickLabel',num2strC(Edg(Tick),'%2.0f'));
xlabel('HTR BAT value');
ylabel('STR BAT value');
FigFN='HTR_STRBATHist';
% saveas(501,[FigFN '.fig']);
% saveas(501,[FigFN '.tif']);
%%
DB=BA-BB-median(BA(FMsk))+median(BB(FMsk));
DB=BA-BB+7;
DB(~FMsk)=0;
Vals=abs(DB(FMsk));
[Ns Bins]=hist(Vals,1000);
figure(505);
plot(Bins,cumsum(Ns)./sum(Ns),'k','LineWidth',2);
set(gca,'FontSize',20);
title(['CDF of BAT difference, median=' num2str(median(Vals))]);
xlabel('BAT difference');
ylabel('%voxels');
FigFN='HRBATCDF';
% saveas(505,[FigFN '.fig']);
% saveas(505,[FigFN '.tif']);
%%
SI=3;
GoodYs=60:200;
GoodXs=50:210;
figure(101010);clf;set(gca,'FontSize',20);
subplot(1,3,1);
imagesc(mritransform(BA(GoodXs,GoodYs,SI)),[0 25]);
set(gca,'FontSize',20);
title('High-temporal resolution');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
colorbar
subplot(1,3,2);
imagesc(mritransform(BB(GoodXs,GoodYs,SI)-7),[0 25]);
set(gca,'FontSize',20);
title('Low-temporal resolution');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
colormap gray
colorbar
subplot(1,3,3);
imagesc(mritransform(abs(DB(GoodXs,GoodYs,SI))),[-5 5])
set(gca,'FontSize',20);
title('Difference');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
colorbar
FigFN='HRBATComparison';
% saveas(101010,[FigFN '.fig']);
% saveas(101010,[FigFN '.tif']);