% 367 FA3D - by fit to B1 from WM=refWMVal
% 368 M03D - Using the same B1, and base M0 from DESPOT1 (Not to be used?)
% 384 DCEM0 - DCEToM0(FinalT1,FA3D);, using baseline, the T1 and FA3D
% 417 BaselineFA3D get FA by the baseline, the FA3D and M03D already calculated
% 429 BaselineFA3Dx get FA by the baseline, the FA3D and DCEM0 already calculated

T12D=NaN([sum(DBrainMask(:)) NumVols]);
for i=1:NumVols
    disp(i);
    CurVol=squeeze(DCE4D(:,:,:,i));
    %     CurT1=DCEVolToT1(FA3D,DCEM0,CurVol);
    CurT1=DCEVolToT1(BaselineFA3Dx,DCEM0,CurVol);
    CurT1(imag(CurT1)~=0)=NaN;
    T12D(:,i)=CurT1(DBrainMask);
end

CTC2DBig=(1./T12D)-repmat(1./FinalT1(DBrainMask),[1 NumVols]);

for qq=SingleAngleIdxs
    dir_to_look_name = dir([WorkingP 'Relaxometry' filesep 'Set_Num_' num2str(qq) filesep 'Core*.nii']);
    TmpFN=[WorkingP 'Relaxometry' filesep 'Set_Num_' num2str(qq) filesep dir_to_look_name(1).name];
    FullVol=loadniidata(TmpFN);
    FullVol4D(:,:,:,qq)=FullVol;
    Tmp=DCEVolToT1(BaselineFA3Dx,DCEM0,FullVol);
    
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
    Tmp=DCEVolToT1(BaselineFA3Dx,DCEM0,FullVol);
    
    FinalT1s(:,:,:,qq)=Tmp;
    Additional_T1_Maps(:,:,:,qq)=Tmp;
    CAdditional_T1_Maps(:,:,:,qq)=Tmp;
end

Additional_T1_Inverse = 1 ./ (FinalT1s(:,:,:,NotBaseIdxs));

% Get the 2D additional maps
AdditionalT1_Inverse_2D = Reshape4d22d(Additional_T1_Inverse,DBrainMask);

Additional_CTC_2D=repPlus(AdditionalT1_Inverse_2D,-R10Clmn);
Additional_CTC_2D(imag(Additional_CTC_2D(:))~=0)=0;
CTC2DBig(imag(CTC2DBig(:))~=0)=0;

Additional_after_main  = Additional_CTC_2D(:,intersect(num_before :end,GoodAnglesF));

%     CTC2DBig = [Additional_before_main CTC2DBig Additional_after_main];
CTC2DBig = [CTC2DBig Additional_after_main];

MskCTCGood=~any(imag(CTC2DBig)~=0,2) | any(isnan(CTC2DBig),2);
CTC2DBigGood=CTC2DBig(MskCTCGood(:),1:UnderSampling:end);

F1=find(DBrainMask);
F2=find(Msk);
B3=ismember(F1,F2);
CTC2DA=CTC2DBig(B3(:),1:UnderSampling:end);

ImagB=any(imag(CTC2DA)~=0,2) | any(isnan(CTC2DA),2);
CTC2DB=CTC2DA(~ImagB,:);
CTC2D=CTC2DB;
if(isempty(CTC2D))
    error('Empty CTC2D');
end

Msk2=Msk;
Msk2(Msk2)=~ImagB;

ImagB3D=Msk;
ImagB3D(Msk)=~ImagB;

[S O]=sort(CTC2D(:,40));
figure;plot(SampleTsNoBefore,CTC2D(O(numel(O)-[4000,5000,6000,100000,50000,3000]),:),'.-')