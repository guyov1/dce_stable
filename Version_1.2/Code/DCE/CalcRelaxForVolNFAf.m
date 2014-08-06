function [ONFAs Mindm MinSdm MinSdmTst]=CalcRelaxForVolNFAf(FA4D,FBrainMaskS,FAs,TRs,nBFGSTries)
%%
FFA2DData=Reshape4d22d(FA4D,FBrainMaskS);
nFAs=numel(FAs);
SFA=size(FA4D);
%%
[AppT1s, AppPDs, AppRMSs]=CalcT1byFAfwSin(FFA2DData',FAs,TRs);
AppT1s3D=FBrainMaskS;
AppT1s3D(FBrainMaskS)=AppT1s<3000;
SAppRMSs=sort(AppRMSs(AppT1s<3000));
T=SAppRMSs(floor(0.995*numel(SAppRMSs)));
AppRMSs3D=FBrainMaskS;
AppRMSs3D(FBrainMaskS)=AppRMSs<T;
FBrainMaskS=AppT1s3D & AppRMSs3D;
% FFA2DData=Reshape4d22d(FA4D,FBrainMaskS);
% nVox=size(FFA2DData,1);
%%
% figure(8723);clf;
MedFA=median(FAs);
[QQQ, MedFAI]=min(abs(FAs-MedFA));
% MedFA=FAs(MedFAI);
NotMedFAI=setdiff(1:nFAs,MedFAI);

nFASmpsF=1000;
% nTestSmps=0;
% nBFGSTries=10;
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
IMsk=FBrainMaskS>100;
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

    IMsk=FBrainMaskS>100;
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
for i=1:nBFGSTries
%     disp(i);
    MinSdmTst(Order(i))=CostFuncnmT(Mindm{Order(i)});
end
[Tmp MinSdmB]=min(MinSdmTst(Order(1:10)));
MinSdmI=Order(MinSdmB);
% BestX0=Mindm{MinSdmI};

    
NFAs=[FAs(NotMedFAI)+Mindm{MinSdmI} FAs(MedFAI)];
ONFAs=FAs;
ONFAs(NotMedFAI)=NFAs(1:numel(NotMedFAI));

disp('End of Detection of change in flip angle');