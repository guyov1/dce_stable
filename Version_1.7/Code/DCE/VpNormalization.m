% WorkingP='C:\STRDCE\John\Database\DCEOut\WhSa_20070813\';
DataP = [WorkingP 'AutoArtBAT' filesep];
D     = dir([WorkingP 'RefVp_*.nii']);
D     = D(strhas({D.name},'WM_'));

if(isempty(D))
    D=dir([WorkingP 'RefAuto_Base_*.nii']);
    D=D(strhas({D.name},'WM_'));
end

RefFN=[WorkingP D(1).name];
AddToLog(WorkingP,'e_t20Norm',['VpNormalization ref file: ' strrep(D(1).name,'_','-')]);

RefVol=loadniidata(RefFN)>0;
Vp3D=loadniidata([DataP 'VpFinal.nii']);
Tmp=Vp3D(RefVol);
Tmp=Tmp(isfinite(Tmp) & Tmp>0);
CurMedWMVp=median(Tmp);
TrgWMVpVal=0.01;
AIFAmpCoeff=CurMedWMVp/TrgWMVpVal;
GeneralDataFN=[WorkingP 'Params.mat'];
save(GeneralDataFN,'AIFAmpCoeff');

%%
RelaxFN=[WorkingP 'Relax.mat'];
load(RelaxFN,'GoodSlices','GoodRows','GoodCols');
WhichModel=loadniidata([DataP 'WhichModel.nii']);

figure(1122);clf;
montage(permute(mritransform(WhichModel(GoodRows,GoodCols,GoodSlices)+1),[1 2 4 3]),[0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0]);
colorbar('YTick',1:5,'YTickLabel',{'No computation','Empty model','Only plasma','Patlak','Tofts'});
set(gcf,'Position',figposition([0 0 100 100]));
title('Model selection');
saveas(1122,[WorkingP 'WhichModel.png']);
saveas(1122,[WorkingP 'WhichModel.fig']);
close(1122);
AddToLog(WorkingP,'ye_t21mdl','WhichModel.','WhichModel.png');
%%
Tmp=(min(1,Vp3D/AIFAmpCoeff));
Tmp(isnan(Vp3D))=NaN;
IRGB3=repmat(permute(Tmp(GoodRows,GoodCols,GoodSlices),[1 2 4 3]),[1 1 1 1]);
figure(1122);clf;
montage(mritransformNoSqueeze(IRGB3*256*10/3),jet(256));
title('{\it{v}}_p range [0 0.3]');
set(gcf,'Position',figposition([0 0 100 100]));

saveas(1122,[WorkingP 'Vp.png']);
saveas(1122,[WorkingP 'Vp.fig']);
close(1122);
AddToLog(WorkingP,'ye_t22vp','{\it{v}}_p','Vp.png');
%%
TmpA=loadniidata([DataP 'KtransFinal.nii']);
Tmp=(min(1,TmpA/AIFAmpCoeff));
Tmp(isnan(TmpA))=NaN;
IRGB3=repmat(permute(Tmp(GoodRows,GoodCols,GoodSlices),[1 2 4 3]),[1 1 1 1]);
figure(1122);clf;
montage(mritransformNoSqueeze(IRGB3*256*10/3),jet(256));
title('{\it{k}}^{trans} range [0 0.3]');
set(gcf,'Position',figposition([0 0 100 100]));

saveas(1122,[WorkingP 'Ktrans.png']);
saveas(1122,[WorkingP 'Ktrans.fig']);
close(1122);
AddToLog(WorkingP,'ye_t23ktrans','{\it{k}}^{trans} range [0 0.3]','Ktrans.png');
%%
TmpA=loadniidata([DataP 'KepFinal.nii']);
Tmp=(min(1,TmpA));
Tmp(isnan(TmpA))=NaN;
IRGB3=repmat(permute(Tmp(GoodRows,GoodCols,GoodSlices),[1 2 4 3]),[1 1 1 1]);
figure(1122);clf;
montage(mritransformNoSqueeze(IRGB3*256),jet(256));
title('{\it{k}}_{ep} range [0 1]');
set(gcf,'Position',figposition([0 0 100 100]));

saveas(1122,[WorkingP 'Kep.png']);
saveas(1122,[WorkingP 'Kep.fig']);
close(1122);
AddToLog(WorkingP,'ye_t24kep','{\it{k}}_{ep} range [0 1]','Kep.png');
%%
Tmp=loadniidata([DataP 'BATFinal.nii']);
Tmp=-Tmp+max(Tmp(:))+1;
IRGB3=repmat(permute(Tmp(GoodRows,GoodCols,GoodSlices),[1 2 4 3]),[1 1 1 1]);
figure(1122);clf;
montage(mritransformNoSqueeze((IRGB3+1)*256/max(Tmp(:))),jet(256));
title('BAT range [0 1]');
set(gcf,'Position',figposition([0 0 100 100]));

saveas(1122,[WorkingP 'BAT.png']);
saveas(1122,[WorkingP 'BAT.fig']);
close(1122);
AddToLog(WorkingP,'ye_t25BAT','BAT range [0 1]','BAT.png');
%%
TmpA=loadniidata([DataP 'VeFinal.nii']);
Tmp=(min(1,TmpA/AIFAmpCoeff));
Tmp(isnan(TmpA))=NaN;
IRGB3=repmat(permute(Tmp(GoodRows,GoodCols,GoodSlices),[1 2 4 3]),[1 1 1 1]);
figure(1122);clf;
montage(mritransformNoSqueeze(IRGB3*256),jet(256));
title('{\it{v}}_{e} range [0 1]');
set(gcf,'Position',figposition([0 0 100 100]));

saveas(1122,[WorkingP 'Ve.png']);
saveas(1122,[WorkingP 'Ve.fig']);
close(1122);
AddToLog(WorkingP,'ye_t26ve','{\it{v}}_{e} range [0 1]','Ve.png');