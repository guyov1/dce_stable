%%
FFA2DData=Reshape4d22d(FA4D,FBrainMaskS);
FFA2DDataF=FFA2DData;
FBrainMaskSF=FBrainMaskS;
tic;[RelaxC3{1:3}]=CalcT1byFAfw2(FFA2DData',FAs,TRs);TRelax(3)=toc;
FRelax4D3=Reshape2DCto4D(RelaxC3,FBrainMaskS);

% FBrainMaskSx=FBrainMaskS>5;
% FBrainMaskSx(:,:,MidSli)=true;
% FFA2DDatax=Reshape4d22d(FA4D,FBrainMaskSx);
% Whichx=[5 6];
% tic;[RelaxC3{1:3}]=CalcT1byFAfw2(FFA2DDatax(:,Whichx)',FAs(Whichx),TRs(Whichx));TRelax(3)=toc;
% FRelax4D3x=Reshape2DCto4D(RelaxC3,FBrainMaskSx);
% figure;imagesc(squeeze(FRelax4D3x(:,:,MidSli,1)),[0 4000]);colorbar;
%%
[AppT1s, AppPDs, AppRMSs]=CalcT1byFAfwSin(FFA2DData',FAs,TRs);
AppT1s3D=FBrainMaskS*0;
AppT1s3D(FBrainMaskS)=AppT1s;
SAppRMSs=sort(AppRMSs(AppT1s<3000));
T=SAppRMSs(floor(0.995*numel(SAppRMSs)));
AppRMSs3D=FBrainMaskS;
AppRMSs3D(FBrainMaskS)=AppRMSs<T;
FBrainMaskS=AppT1s3D<3000 & AppRMSs3D;
FFA2DData=Reshape4d22d(FA4D,FBrainMaskS);
nVox=size(FFA2DData,1);
%%
MaxT1Disp=2500;
MaxRMSDisp=40;
% figure(8723);clf;
% subplot(1,2,1);imagesc(squeeze(mritransform(FRelax4D3(:,:,CurSliI,1))),[0 MaxT1Disp]);
% subplot(1,2,2);imagesc(squeeze(mritransform(FRelax4D3(:,:,CurSliI,3))),[0
% MaxRMSDisp]);title(num2str(RRMS(3)))
MedFA=median(FAs);
[QQQ, MedFAI]=min(abs(FAs-MedFA));
MedFA=FAs(MedFAI);
NotMedFAI=setdiff(1:nFAs,MedFAI);

nFASmpsF=1000;
nTestSmps=0;
nBFGSTries=10;
Mindm=cell(1,nBFGSTries);
MinSdm=zeros(1,nBFGSTries);
MinSdmTst=zeros(1,nBFGSTries);
xlFAd=-ones(1,nFAs-1).*min(11,FAs(NotMedFAI)-0.5);
xhFAd=ones(1,nFAs-1).*11;
x0FAd=zeros(1,nFAs-1);
MaskForSamples=FBrainMaskS;
MaskForSamples(:,:,1:(floor(SFA(3)/2)-2))=false;
MaskForSamples(:,:,(ceil(SFA(3)/2)+2):end)=false;
% Idxs=getKrandomSamples(find(MaskForSamples),nFASmpsF*nBFGSTries + nTestSmps);
Idxs=find(MaskForSamples);
% Idx=find(MaskForSamples);
IMsk=FBrainMask>100;
IMsk(Idxs)=true;
FAData=Reshape4d22d(FA4D,IMsk)';
MFADataT=[FAData(NotMedFAI,:); FAData(MedFAI,:)];
RelaxFMT=@(x) CalcT1byFAfw2(MFADataT,x,TRs);
CostFuncnmT=@(x) mean(getKthOutput(3,RelaxFMT,{[FAs(NotMedFAI) FAs(MedFAI)]+[x 0]}));

disp(['Start of Detection of change in flip angle, nFASmpsF = ' num2str(nFASmpsF)]);
%%
for i=1:nBFGSTries
    disp(i);
%     Idx=Idxs(nTestSmps+nFASmpsF*(i-1)+(1:nFASmpsF));
    Idx=getKrandomSamples(Idxs,nFASmpsF);

    IMsk=FBrainMask>100;
    IMsk(Idx)=true;
    FAData=Reshape4d22d(FA4D,IMsk)';
        
    MFAData=[FAData(NotMedFAI,:); FAData(MedFAI,:)];
    RelaxFM=@(x) CalcT1byFAfw2(MFAData,x,TRs);
    CostFuncnm=@(x) mean(getKthOutput(3,RelaxFM,{[FAs(NotMedFAI) FAs(MedFAI)]+[x 0]}));
    % Detection of change in flip angle
    Mindm{i}=BFGSf(CostFuncnm,[x0FAd ],[xlFAd ],[xhFAd ],struct('GradObj','off','Display','off'));
    MinSdm(i)=CostFuncnm(Mindm{i});
%     MinSdmTst(i)=CostFuncnmT(Mindm{i});
end
[BestMin1 Order]=sort(MinSdm);
for i=1:10
%     disp(i);
    MinSdmTst(Order(i))=CostFuncnmT(Mindm{Order(i)});
end
[Tmp MinSdmB]=min(MinSdmTst(Order(1:10)));
MinSdmI=Order(MinSdmB);
BestX0=Mindm{MinSdmI};

% for i=1:nBFGSTries
%     disp(i);
% %     Idx=Idxs(nTestSmps+nFASmpsF*(i-1)+(1:nFASmpsF));
%     Idx=getKrandomSamples(Idxs,nFASmpsF);
% 
%     IMsk=FBrainMask>100;
%     IMsk(Idx)=true;
%     FAData=Reshape4d22d(FA4D,IMsk)';
%         
%     MFAData=[FAData(NotMedFAI,:); FAData(MedFAI,:)];
%     RelaxFM=@(x) CalcT1byFAfw2(MFAData,x,TRs);
%     CostFuncnm=@(x) mean(getKthOutput(3,RelaxFM,{[FAs(NotMedFAI) FAs(MedFAI)]+[x 0]}));
%     % Detection of change in flip angle
%     Mindm{i}=BFGSf(CostFuncnm,[BestX0 ],[xlFAd ],[xhFAd ],struct('GradObj','off','Display','off'));
%     MinSdm(i)=CostFuncnm(Mindm{i});
% %     MinSdmTst(i)=CostFuncnmT(Mindm{i});
% end
% [BestMin1 Order]=sort(MinSdm);
% for i=1:10
%     disp(i);
%     MinSdmTst(Order(i))=CostFuncnmT(Mindm{Order(i)});
% end
% [Tmp MinSdmB]=min(MinSdmTst(Order(1:10)));
% MinSdmI=Order(MinSdmB);
% Mindm{MinSdmI}
% gDispMat([MinSdm;MinSdmTst]);

% A=cat(1,Mindm{:});
% B=zeros(nBFGSTries,nFAs);
% B(:,NotMedFAI)=repmat(FAs(NotMedFAI),nBFGSTries,1) + A;
% B(:,MedFAI)=MedFA;
% CurRFAs=RFAs(TRGroupC{CurCGroup});
% NCurRFAs=CurRFAs*MedFA/CurRFAs(MedFAI);
% BestD=NCurRFAs(NotMedFAI)-FAs(NotMedFAI);
% C=B*CurRFAs(MedFAI)/MedFA;
% BestCostA=CostFuncnm(BestD);
% BestCostF=CostFuncnmT(BestD);
% gDispMat([[BestCostA MinSdm];[BestCostF MinSdmTst]]);
% C(Order(1:10),:)
% C(MinSdmI,:)
%%
% [Tmp MinSdmI]=min(MinSdm);
% disp(MinSdm);
% disp(['NFAs Found by ' num2str(MinSdmI) ', min= ' num2str(Tmp)]);
    
NFAs=[FAs(NotMedFAI)+Mindm{MinSdmI} FAs(MedFAI)];
ONFAs=FAs;
ONFAs(NotMedFAI)=NFAs(1:numel(NotMedFAI));

disp('End of Detection of change in flip angle');