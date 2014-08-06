% WorkingP='C:\STRDCE\John\Database\DCEOut\WhSa_20070813\';
% DataP = [WorkingP 'AutoArtBAT' filesep];
D     = dir([WorkingP 'RefVp_*.nii']);
D     = D(strhas({D.name},'WM_'));

if(isempty(D))
    D=dir([WorkingP 'RefAuto_Base_*.nii']);
    D=D(strhas({D.name},'WM_'));
end

RefFN=[WorkingP D(1).name];
AddToLog(WorkingP,['e' VPstr '_t20Norm'],[VPstr ' VpNormalization ref file: ' strrep(D(1).name,'_','-')]);

RefVol=loadniidata(RefFN)>0;
Vp3D=loadniidata([DataP 'VpFinal.nii']);
Tmp=Vp3D(RefVol);
Tmp=Tmp(isfinite(Tmp) & Tmp>0);
CurMedWMVp=median(Tmp);
TrgWMVpVal=0.01;
AIFAmpCoeff=CurMedWMVp/TrgWMVpVal;
GeneralDataFN=[WorkingP 'Params.mat'];
save(GeneralDataFN,'AIFAmpCoeff');

MeanFN=[WorkingP 'DCEMean.nii'];
Raw2Nii(Vp3D/AIFAmpCoeff,[DataP 'VpFinalN.nii'],'float32', MeanFN);

TmpA=loadniidata([DataP 'Vp1.nii']);
Raw2Nii(TmpA/AIFAmpCoeff,[DataP 'Vp1N.nii'],'float32', MeanFN);
TmpA=loadniidata([DataP 'Vp2.nii']);
Raw2Nii(TmpA/AIFAmpCoeff,[DataP 'Vp2N.nii'],'float32', MeanFN);
TmpA=loadniidata([DataP 'Vp3.nii']);
Raw2Nii(TmpA/AIFAmpCoeff,[DataP 'Vp3N.nii'],'float32', MeanFN);
%%
RelaxFN=[WorkingP 'Relax.mat'];
load(RelaxFN,'GoodSlices','GoodRows','GoodCols');
WhichModel=loadniidata([DataP 'WhichModel.nii']);

figure(1122);clf;
montage(permute(mritransform(WhichModel(GoodRows,GoodCols,GoodSlices)+1),[1 2 4 3]),[0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0]);
colorbar('YTick',1:5,'YTickLabel',{'No computation','Empty model','Only plasma','Patlak','Tofts'});
set(gcf,'Position',figposition([0 0 100 100]));
title('Model selection');
saveas(1122,[WorkingP VPstr 'WhichModel.png']);
saveas(1122,[WorkingP VPstr 'WhichModel.fig']);
close(1122);
AddToLog(WorkingP,['ye' VPstr '_t21mdl'],[VPstr ' WhichModel.'],[VPstr 'WhichModel.png']);
%%
Tmp=(min(1,Vp3D/AIFAmpCoeff));
Tmp(isnan(Vp3D))=NaN;
IRGB3=repmat(permute(Tmp(GoodRows,GoodCols,GoodSlices),[1 2 4 3]),[1 1 1 1]);
figure(1122);clf;
montage(mritransformNoSqueeze(IRGB3*256*10/3),jet(256));
title('{\it{v}}_p range [0 0.3]');
set(gcf,'Position',figposition([0 0 100 100]));

saveas(1122,[WorkingP VPstr 'Vp.png']);
saveas(1122,[WorkingP VPstr 'Vp.fig']);
close(1122);
AddToLog(WorkingP,['ye' VPstr '_t22vp'],[VPstr ' {\it{v}}_p'],[VPstr 'Vp.png']);
%%
TmpA=loadniidata([DataP 'KtransFinal.nii']);
Tmp=(min(1,TmpA/AIFAmpCoeff));
Tmp(isnan(TmpA))=NaN;
IRGB3=repmat(permute(Tmp(GoodRows,GoodCols,GoodSlices),[1 2 4 3]),[1 1 1 1]);
figure(1122);clf;
montage(mritransformNoSqueeze(IRGB3*256*10/3),jet(256));
title('{\it{k}}^{trans} range [0 0.3]');
set(gcf,'Position',figposition([0 0 100 100]));

saveas(1122,[WorkingP VPstr 'Ktrans.png']);
saveas(1122,[WorkingP VPstr 'Ktrans.fig']);
close(1122);
AddToLog(WorkingP,['ye' VPstr '_t23ktrans'],[VPstr ' {\it{k}}^{trans} range [0 0.3]'],[VPstr 'Ktrans.png']);

Raw2Nii(TmpA/AIFAmpCoeff,[DataP 'KtransFinalN.nii'],'float32', MeanFN);

TmpA=loadniidata([DataP 'Ktrans3.nii']);
Raw2Nii(TmpA/AIFAmpCoeff,[DataP 'Ktrans3N.nii'],'float32', MeanFN);
TmpA=loadniidata([DataP 'Ktrans2.nii']);
Raw2Nii(TmpA/AIFAmpCoeff,[DataP 'Ktrans2N.nii'],'float32', MeanFN);
%%
TmpA=loadniidata([DataP 'Kepfinal.nii']);
Tmp=(min(1,TmpA));
Tmp(isnan(TmpA))=NaN;
IRGB3=repmat(permute(Tmp(GoodRows,GoodCols,GoodSlices),[1 2 4 3]),[1 1 1 1]);
figure(1122);clf;
montage(mritransformNoSqueeze(IRGB3*256),jet(256));
title('{\it{k}}_{ep} range [0 1]');
set(gcf,'Position',figposition([0 0 100 100]));

saveas(1122,[WorkingP VPstr 'Kep.png']);
saveas(1122,[WorkingP VPstr 'Kep.fig']);
close(1122);
AddToLog(WorkingP,['ye' VPstr '_t24kep'],[VPstr ' {\it{k}}_{ep} range [0 1]'],[VPstr 'Kep.png']);
%%
Tmp=loadniidata([DataP 'BATfinal.nii']);
Tmp=-Tmp+max(Tmp(:))+1;
IRGB3=repmat(permute(Tmp(GoodRows,GoodCols,GoodSlices),[1 2 4 3]),[1 1 1 1]);
figure(1122);clf;
montage(mritransformNoSqueeze((IRGB3+1)*256/max(Tmp(:))),jet(256));
title('BAT range [0 1]');
set(gcf,'Position',figposition([0 0 100 100]));

saveas(1122,[WorkingP VPstr 'BAT.png']);
saveas(1122,[WorkingP VPstr 'BAT.fig']);
close(1122);
AddToLog(WorkingP,['ye' VPstr '_t25BAT'],[VPstr ' BAT range [0 1]'],[VPstr 'BAT.png']);
%%
TmpA=loadniidata([DataP 'VeFinal.nii']);
Tmp=(min(1,TmpA/AIFAmpCoeff));
Tmp(isnan(TmpA))=NaN;
IRGB3=repmat(permute(Tmp(GoodRows,GoodCols,GoodSlices),[1 2 4 3]),[1 1 1 1]);
figure(1122);clf;
montage(mritransformNoSqueeze(IRGB3*256),jet(256));
title('{\it{v}}_{e} range [0 1]');
set(gcf,'Position',figposition([0 0 100 100]));

saveas(1122,[WorkingP VPstr 'Ve.png']);
saveas(1122,[WorkingP VPstr 'Ve.fig']);
close(1122);
AddToLog(WorkingP,['ye' VPstr '_t26ve'],[VPstr ' {\it{v}}_{e} range [0 1]'],[VPstr 'Ve.png']);