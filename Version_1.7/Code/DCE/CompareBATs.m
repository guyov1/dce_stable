FA=loadniidata('C:\DCE\John\Database\DCEOut\SmVl_20120930\F.nii');
FB=loadniidata('C:\DCE\John\Database\DCEOut\SmVl_20120930\F_3.nii');
RMS=loadniidata('C:\DCE\John\Database\DCEOut\SmVl_20120930\RMS.nii');
FThresh=10^-3;
FMsk=FA<FThresh & FB<FThresh & RMS>0;
figure;imagesc(mritransform(FMsk(:,:,2)));
%%
BA=loadniidata('C:\DCE\John\Database\DCEOut\SmVl_20120930\BAT.nii');
BB=loadniidata('C:\DCE\John\Database\DCEOut\SmVl_20120930\BAT_3.nii');

Edg=unique(BA(FMsk));
HH=hist3([BA(FMsk) BB(FMsk)],{Edg,Edg});
figure(501);clf;
set(gca,'FontSize',20);
axis equal
imagesc(HH);
[R P]=corrcoef([BA(FMsk) BB(FMsk)]);
% title(getKthElement(corrcoef([BA(FMsk) BB(FMsk)]),2));
title(['HTR vs. STR 2D histogram, correlation=' num2str(R(2),'%2.2g'), ', p<1e-10']);
colormap gray
Tick=floor(linspace(1,ceil(numel(Edg)/2),3));
Tick=[Tick numel(Edg)-Tick(2)+1 numel(Edg)];
set(gca,'XTick',Tick);
set(gca,'YTick',Tick);
set(gca,'XTickLabel',num2strC(Edg(Tick),'%2.1g'));
set(gca,'YTickLabel',num2strC(Edg(Tick),'%2.1g'));
xlabel('HTR BAT value');
ylabel('STR BAT value');
FigFN='HTR_STRBATHist';
saveas(501,[FigFN '.fig']);
saveas(501,[FigFN '.tif']);
%%
DB=BA-BB-median(BA(FMsk))+median(BB(FMsk));
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
saveas(505,[FigFN '.fig']);
saveas(505,[FigFN '.tif']);
%%
SI=3;
GoodYs=60:200;
GoodXs=50:210;
figure(101010);clf;set(gca,'FontSize',20);
subplot(1,3,1);
imagesc(mritransform(BA(GoodXs,GoodYs,SI)));
set(gca,'FontSize',20);
title('High-temporal resolution');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
colorbar
subplot(1,3,2);
imagesc(mritransform(BB(GoodXs,GoodYs,SI)));
set(gca,'FontSize',20);
title('Low-temporal resolution');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
colormap gray
colorbar
subplot(1,3,3);
imagesc(mritransform(abs(DB(GoodXs,GoodYs,SI))),[-1 1])
set(gca,'FontSize',20);
title('Difference');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
colorbar
FigFN='HRBATComparison';
saveas(101010,[FigFN '.fig']);
saveas(101010,[FigFN '.tif']);