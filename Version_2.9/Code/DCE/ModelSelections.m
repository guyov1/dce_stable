% WorkingP='C:\STRDCE\John\Database\DCEOut\RuYe_20091224\';
% DataP=[WorkingP 'AutoArtBAT' filesep];
% a=load([WorkingP 'PKM3D.mat']);
RSSIdxs=7:10;
NParams=[0 2 3 4];
%%
% MAD=median(abs(diff(CTC2DBigGood,[],2)),2);
% SSig=(MAD/0.67).^2;.*repmat(SSig,[1 4])
AIC=(0*log(nVols)+1*2)*repmat(NParams,[size(PKs,1) 1])+nVols*0.1*log(PKs(:,RSSIdxs));
Correction=2*NParams.*(NParams+1)./(nVols-NParams-1);
AICc=AIC+repmat(Correction,[size(PKs,1) 1]);

[Tmp, ChosenByAIC]=min(AICc,[],2);

Tmp3D(MskX)=ChosenByAIC;
figure(1122);clf;
montage(permute(mritransform(Tmp3D(GoodRows,GoodCols,GoodSlices)+1),[1 2 4 3]),[0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0]);
colorbar('YTick',1:5,'YTickLabel',{'No computation','Empty model','Only plasma','Patlak','Tofts'});
set(gcf,'Position',figposition([0 0 100 100]));
title('AICc');
saveas(1122,[WorkingP VPstr 'AICc.png']);
saveas(1122,[WorkingP VPstr 'AICc.fig']);
close(1122);
AddToLog(WorkingP,['ye' VPstr '_t21AICc'],[VPstr ' AICc.'],[VPstr 'AICc.png']);
%%
Raw2Nii(Tmp3D,[PKOutP 'AICc' '.nii'],'float32', MeanFN);