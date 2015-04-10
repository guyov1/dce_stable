for qq=SingleAngleIdxs
    dir_to_look_name = dir([WorkingP 'Relaxometry' filesep 'Set_Num_' num2str(qq) filesep 'Core*.nii']);
    TmpFN=[WorkingP 'Relaxometry' filesep 'Set_Num_' num2str(qq) filesep dir_to_look_name(1).name];
    FullVol=loadniidata(TmpFN);
    FullVol4D(:,:,:,qq)=FullVol;
    Tmp=DCEVolToT1(FA3D,DCEM0,FullVol);
    
    FinalT1s(:,:,:,qq)=Tmp;
    Additional_T1_Maps(:,:,:,qq)=Tmp;
    CAdditional_T1_Maps(:,:,:,qq)=Tmp;
end

for qq=DESPOT1Idxs
    dir_to_look_name = dir([WorkingP 'Relaxometry' filesep 'Set_Num_' num2str(qq) filesep 'Core*.nii']);
    Tmp=cell2mat({dir_to_look_name.name}');
    TmpFAs=str2num(Tmp(:,16:17));
    [AA,Chosen]=min(abs(TmpFAs-DCEFA));
    TmpFN=[WorkingP 'Relaxometry' filesep 'Set_Num_' num2str(qq) filesep dir_to_look_name(Chosen).name];
    FullVol=loadniidata(TmpFN);
    FullVol4D(:,:,:,qq)=FullVol;
    Tmp=DCEVolToT1(FA3D,DCEM0,FullVol);
    
    FinalT1s(:,:,:,qq)=Tmp;
    Additional_T1_Maps(:,:,:,qq)=Tmp;
    CAdditional_T1_Maps(:,:,:,qq)=Tmp;
end

[S O]=sort(CTC2D(:,40));
%%
T1Mx=3000;
figure(232);clf;
subplot(2,2,1);
Tmp(ImagB3D)=CTC2D(:,NumVols);
Tmp=1./(Tmp+R10);
Tmp(imag(Tmp)~=0)=0;
imagesc(mritransform(Tmp(:,:,MidSli)),[0 T1Mx])
title('Last of DCEMain');
subplot(2,2,2);
Tmp=DCEVolToT1(FA3D,DCEM0,squeeze(FullVol4D(:,:,:,SingleAngleIdxs(end))));
Tmp(imag(Tmp)~=0)=0;
TmpB=Tmp;
imagesc(mritransform(Tmp(:,:,MidSli)),[0 T1Mx])
title('Last Single');
subplot(2,2,3);
Tmp=DCEVolToT1(FA3D,DCEM0,squeeze(FullVol4D(:,:,:,DESPOT1Idxs(end))));
Tmp(imag(Tmp)~=0)=0;
TmpC=Tmp;
imagesc(mritransform(Tmp(:,:,MidSli)),[0 T1Mx])
title('Last DESPOT1');
subplot(2,2,4);
Tmp=TmpB-TmpC;
imagesc(mritransform(Tmp(:,:,MidSli)),[-100 100]);