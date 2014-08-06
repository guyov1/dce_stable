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

nFASmpsF=100;
nTestSmps=0;
nBFGSTries=3;
Mindm=cell(1,nBFGSTries);
MinSdm=zeros(1,nBFGSTries);
MinSdmTst=zeros(1,nBFGSTries);
xlFAd=-ones(1,nFAs-1).*min(11,FAs(NotMedFAI)-0.2);
xhFAd=ones(1,nFAs-1).*21;
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
%
MindmM=cell(1,nBFGSTries);
MinSdmM=zeros(1,nBFGSTries);
MinSdmTstM=zeros(1,nBFGSTries);

FAsX=[FAs(NotMedFAI) FAs(MedFAI)];
TRsX=[TRs(NotMedFAI) TRs(MedFAI)];
for i=1:nBFGSTries
    disp(i);
%     Idx=Idxs(nTestSmps+nFASmpsF*(i-1)+(1:nFASmpsF));
    Idx=getKrandomSamples(Idxs,nFASmpsF);

    IMsk=FBrainMask>100;
    IMsk(Idx)=true;
    FAData=Reshape4d22d(FA4D,IMsk)';
        
    MFAData=[FAData(NotMedFAI,:); FAData(MedFAI,:)];
    RelaxFM=@(x) CalcT1byFAfw2(MFAData,x,TRsX);
    CostFuncnm=@(x) mean(getKthOutput(3,RelaxFM,{[FAs(NotMedFAI) FAs(MedFAI)]+[x 0]}));
    
    RelaxFMMul=@(x,y) CalcT1byFAfw2(MFAData.*repmat(([1 y]/min([1 y]))',1,size(MFAData,2)),x,TRs);
    CostFuncnmMul=@(x) mean(getKthOutput(3,RelaxFMMul,{[FAs(NotMedFAI) FAs(MedFAI)]+[x(1:(nFAs-1)) 0],x(nFAs:end)}));
    
    % Detection of change in flip angle
    Mindm{i}=BFGSf(CostFuncnm,[x0FAd ],[xlFAd ],[xhFAd ],struct('GradObj','off','Display','off'));
    MinSdm(i)=CostFuncnm(Mindm{i});
%     MinSdmTst(i)=CostFuncnmT(Mindm{i});
    disp('b');
    MindmM{i}=BFGSf(CostFuncnmMul,[x0FAd x0FAd+1],[xlFAd x0FAd+0.5],[xhFAd x0FAd+2],struct('GradObj','off','Display','off'));
    MinSdmM(i)=CostFuncnmMul(MindmM{i});

%     x=x0FAd;
%     for tt=1:10
%         disp(tt);
%         [Tmp{1:3}]=CalcT1byFAfw2(MFAData,[FAs(NotMedFAI) FAs(MedFAI)]+[x 0],TRsX);
%         for jj=1:nFAs-1
%             CurCost=@(y) gCost(MFAData(jj,:)',SPGRfM(Tmp{1}',Tmp{2}',y,TRsX(jj)),'RMS');
%             x(jj)=BFGSf(CurCost,FAsX(jj)+x(jj),0.5,FAsX(jj)+11,struct('GradObj','off','Display','off'))-FAsX(jj);
%         end
%         disp(CostFuncnm(x));
%     end
%     SimF=SPGRfM(Tmp{1}',Tmp{2}',[FAs(NotMedFAI) FAs(MedFAI)]+[x 0],TRs);
end
%%
[BestMin1 Order]=sort(MinSdm);
for i=1:nBFGSTries
%     disp(i);
    MinSdmTst(Order(i))=CostFuncnmT(Mindm{Order(i)});
end
[Tmp MinSdmB]=min(MinSdmTst(Order(1:nBFGSTries)));
MinSdmI=Order(MinSdmB);
BestX0=Mindm{MinSdmI};

    
NFAs=[FAs(NotMedFAI)+Mindm{MinSdmI} FAs(MedFAI)];
ONFAs=FAs;
ONFAs(NotMedFAI)=NFAs(1:numel(NotMedFAI));

disp('End of Detection of change in flip angle');
%%
FFA2DData=Reshape4d22d(FA4D(:,:,:,:),FBrainMaskS);
[RelaxCG{1:3}]=CalcT1byFAfw2(FFA2DData(:,GoodFAsF)',ONFAs(GoodFAsF),TRs(GoodFAsF));

CostFunc=@(FA) gCost(FFA2DData(:,end),SPGRfM(RelaxCG{1}',RelaxCG{2}',FA,TRs(end)),'RMS');
BestBaseFA=BFGSf(CostFunc,FAs(end),0.5,FAs(end)*3,struct('GradObj','off','Display','off'));

SimByNomFA=SPGRfM(RelaxCG{1}',RelaxCG{2}',FAs(end),TRs(end));
CostFunc=@(k) gCost(FFA2DData(:,end),SimByNomFA*k,{'RMS'});
BestBaseFac=BFGSf(CostFunc,1,0.01,10,struct('GradObj','off','Display','off'));

%%
Data=FBrainMaskS*0;
Data(FBrainMaskS)=FFA2DData(:,end);

SimByNominalFA=FBrainMaskS*0;
SimByNominalFA(FBrainMaskS)=SPGRfM(RelaxCG{1}',RelaxCG{2}',FAs(end),TRs(end));

SimByBestBase=FBrainMaskS*0;
SimByBestBase(FBrainMaskS)=SPGRfM(RelaxCG{1}',RelaxCG{2}',BestBaseFA,TRs(end));


[Mn Mx]=FindDR(FFA2DData(:,end));
figure(191919);clf;
subplot(2,3,1);
imagesc(mritransform(Data(:,:,MidSli)),[0 Mx]);
subplot(2,3,2);
imagesc(mritransform(SimByNominalFA(:,:,MidSli)*BestBaseFac),[0 Mx]);
subplot(2,3,3);
imagesc(mritransform(SimByBestBase(:,:,MidSli)),[0 Mx]);
subplot(2,3,5);
imagesc(mritransform(Data(:,:,MidSli))-mritransform(SimByNominalFA(:,:,MidSli)*BestBaseFac),[-Mx/10 Mx/10]);
subplot(2,3,6);
imagesc(mritransform(Data(:,:,MidSli))-mritransform(SimByBestBase(:,:,MidSli)),[-Mx/10 Mx/10]);
%%
FFA2DData=Reshape4d22d(FA4D(:,:,:,:),FBrainMaskS);
[RelaxCG{1:3}]=CalcT1byFAfw2(FFA2DData(:,GoodFAsF)',ONFAs(GoodFAsF),TRs(GoodFAsF));
%%
Tmp=FBrainMaskS*0;
figure(8181);clf;
for i=1:nFAs
    gsubplot(3,nFAs,1,i);
    Tmp(FBrainMaskS)=FFA2DData(:,i);
    [Mn Mx]=FindDR(FFA2DData(:,i));
    imagesc(mritransform(Tmp(:,:,MidSli)),[0 Mx]);
    gsubplot(3,nFAs,2,i);
    Tmp(FBrainMaskS)=SPGRfM(RelaxCG{1}',RelaxCG{2}',ONFAs(i),TRs(end));
    imagesc(mritransform(Tmp(:,:,MidSli)),[0 Mx]);
    gsubplot(3,nFAs,3,i);
    Tmp(FBrainMaskS)=FFA2DData(:,i)-Tmp(FBrainMaskS);
    imagesc(mritransform(Tmp(:,:,MidSli)),[-100 100]);
end