s = RandStream('mt19937ar','Seed',1837463.8374526);
RandStream.setGlobalStream(s);
warning('off','MATLAB:rankDeficientMatrix');

N=1000;
nIter=250;
NoiseSig=0.05;
FANominal=[3 10 15 20 30];
FATrue=[2.59 11.31 15 22.29 28.53];
% FATrue=FANominal;
% Gains=[1 1 1 1 1];
Gains=[1.27 0.95 1 0.83 2];
% TRs=[5.68 5.68 5.33 5.11 5.02];
TRs=[5.68 5.68 5.68 5.68 5.68];

M0=10000*(1+randn(1,N)*0.3);
T1=500+rand(1,N)*1500;
nFAs=numel(FANominal);
SNoNoise=SPGRfM(T1,M0,FATrue,TRs);
SNoNoiseGain=repMulti(SNoNoise,Gains);
S=SNoNoiseGain.*(1+randn(N,nFAs)*NoiseSig);
% Usual T1 Estimation
[RelaxBase{1:3}]=CalcT1byFAfw2(S',FANominal,TRs);
% Init
CurGains=ones(1,nFAs);
CurFAs=FANominal;
CurRelax=RelaxBase;
CurTRFac=S*0+1;
MedNomFA=median(FANominal);
MedNomFAI=find(FANominal==MedNomFA,1);
disp('--');
Iter=1;
%%
Others=setdiff(1:nFAs,[1 MedNomFAI 5]);
FAOrder=[1 MedNomFAI 5 Others];
% T1M0From2(S(:,FAOrder), FANominal(FAOrder)+[0 0 zeros(1,nFAs-2)], TRs(FAOrder));
CostFuncByTwo=@(x) T1M0From4(repMulti(S(:,FAOrder),[1/x(3) 1 1/x(4) ones(1,nFAs-3)]), FANominal(FAOrder)+[x(1) 0 x(2) zeros(1,nFAs-3)], TRs(FAOrder));

x0=[0 0 1 1];
xReal=[FATrue([1 5])-FANominal([1 5]) Gains([1 5])];

CostFuncByTwo(x0)
CostFuncByTwo(xReal)

disp('-1-');
xBest=fminsearch(CostFuncByTwo,x0);
disp('-2-');
[CostFuncByTwo(x0) CostFuncByTwo(xBest) CostFuncByTwo(xReal)]
[x0; xBest; xReal]
%%
Others=setdiff(1:nFAs,[1 MedNomFAI]);
FAOrder=[1 MedNomFAI Others];
% T1M0From2(S(:,FAOrder), FANominal(FAOrder)+[0 0 zeros(1,nFAs-2)], TRs(FAOrder));
CostFuncByOne=@(x) T1M0From2(repMulti(S(:,FAOrder),[1/x(2) ones(1,nFAs-1)]), FANominal(FAOrder)+[x(1) zeros(1,nFAs-1)], TRs(FAOrder));

x0=[0 1];
xReal=[FATrue(1)-FANominal(1) Gains(1)];

CostFuncByOne(x0)
CostFuncByOne(xReal)

disp('-1-');
xBest=fminsearch(CostFuncByOne,x0);
disp('-2-');
[CostFuncByOne(x0) CostFuncByOne(xBest) CostFuncByOne(xReal)]
[x0; xBest; xReal]
%%
% % x0=[0 0 0 0];
% WhichFAs=1:5;
% x0=[CurFAs([1 2 4 5])-FANominal([1 2 4 5]) CurGains([1 2 4 5])-[1 1 1 1]];
% CostFunc=@(x) getKthElement(sort(getKthOutput(3,@CalcT1byFAfw2,{repMulti(S(:,WhichFAs),1./[min(3,max(0.3,1+x(5))) min(3,max(0.3,1+x(6))) 1 min(3,max(0.3,1+x(7))) min(3,max(0.3,1+x(8)))])',FANominal(WhichFAs)+[x(1) x(2) 0 x(3) x(4)],TRs(WhichFAs)})),floor(N*0.5));
% % CostFunc2=@(x) mean(getKthOutput(3,@CalcT1byFAfw2,{repMulti(S(:,[1 3 5]),1./[min(3,max(0.3,1+x(3))) 1 min(3,max(0.3,1+x(4)))])',FANominal([1 3 5])+[x(1) 0 x(2)],TRs([1 3 5])}));
% disp('a');
% BestX=fminsearch(CostFunc,x0)
% x=BestX;
% xReal=[FATrue([1 2 4 5])-FANominal([1 2 4 5]) Gains([1 2 4 5])-[1 1 1 1]]
% [CostFunc(BestX) CostFunc(xReal) CostFunc(x0)]
% % [CostFunc2(BestX) CostFunc2(xReal) CostFunc2(x0)]
% [x0; BestX; xReal]
%%
for Iter=1:nIter
    disp(Iter);
%     for ii=1:nFAs
%         FARange=max(0.5,FANominal(ii)-10):0.1:(FANominal(ii)+10);
%         FAsSim=SPGRfM(CurRelax{1}',CurRelax{2}',FARange,FARange*0+TRs(ii));
%         clear GainForii
%         for j=1:numel(FARange)
%             GainForii(j)=FAsSim(:,j)\S(:,ii);
%         end
%         [MinSRMS,MIdx]=min(mean(repPlus(repMulti(FAsSim',GainForii'),-S(:,ii)').^2,2));
%         CurFAs(ii)=FARange(MIdx);
%         CurGains(ii)=GainForii(MIdx);
%     end
%     GainFactor=1/CurGains(MedNomFAI);
%     CurGains=CurGains*GainFactor;
%     FAFactor=MedNomFA/CurFAs(MedNomFAI);
%     CurFAs=CurFAs*FAFactor;
%     AllFAs(Iter,:)=CurFAs;
%     AllGains(Iter,:)=CurGains;
% if(1)
%     Estimate Gains by current state
    CurSim=SPGRfM(CurRelax{1}',CurRelax{2}',CurFAs,TRs);
    OldGains=CurGains;
    CurGains=median(S./CurSim,1);
    GainFactor=1/CurGains(MedNomFAI);
    CurGains=CurGains*GainFactor;
    AllGains(Iter,:)=CurGains;
%     [CurRelax{1:3}]=CalcT1byFAfw2(repMulti(S,1./CurGains)',CurFAs,TRs);
%     Estimate FA by current state
    for ii=1:nFAs
        CurFAs(ii)=CalcFAgSigT1M0(S(:,ii)/CurGains(ii),CurRelax{1},CurRelax{2},TRs(ii),100,max(0.5,FANominal(ii)/2),FANominal(ii)*2);
    end
    FAFactor=MedNomFA/CurFAs(MedNomFAI);
    CurFAs=CurFAs*FAFactor;
    AllFAs(Iter,:)=CurFAs;
%     Find TR compensation factors
    CurSim=SPGRfM(CurRelax{1}',CurRelax{2}',CurFAs,TRs);
    SimTR1=SPGRfM(CurRelax{1}',CurRelax{2}',CurFAs,ones(1,nFAs)*TRs(1));
    CurTRFac=CurSim./SimTR1;
%     Reestimate M0,T1
    [CurRelax{1:3}]=CalcT1byFAfw2(repMulti(S./CurTRFac,1./CurGains)',CurFAs,TRs);
    AllRelax{Iter}=CurRelax;
    rT1ErrorsCur=log2(CurRelax{1}'./T1);
    Xs=-1:0.01:1;
    NECur=histc(rT1ErrorsCur,Xs);
    T1acc(Iter)=2.^Xs(find(NECur==max(NECur),1));
end
%%
for Iter=1:nIter
    TT=sort(AllRelax{Iter}{3});
    RMS(Iter)=TT(floor(numel(TT)*0.8));
%     RMS(Iter)=median(AllRelax{Iter}{3});
end
[MinRMS MinRMSI]=min(RMS);
CurFAs=AllFAs(MinRMSI,:);
CurGains=AllGains(MinRMSI,:);
CurRelax=AllRelax{MinRMSI};

[Gains; CurGains; FATrue; CurFAs]
rT1ErrorsBase=log2(RelaxBase{1}'./T1);
rT1ErrorsCur=log2(CurRelax{1}'./T1);
Xs=-1:0.01:1;
NEBase=histc(rT1ErrorsBase,Xs);
NECur=histc(rT1ErrorsCur,Xs);
figure(828384);clf;
subplot(2,2,1);
plot(Xs,NEBase./sum(NEBase),'r',Xs,NECur./sum(NECur),'b');
Tmp=get(gca,'XTickLabel');
set(gca,'XTickLabel',num2str(2.^str2num(Tmp),'%2.2f'));
title([2.^Xs(find(NEBase==max(NEBase),1)) 2.^Xs(find(NECur==max(NECur),1))]);

subplot(2,2,2);
plot(log2(repMulti(AllFAs,1./FATrue)));hold on;
plot(log2(repMulti(AllGains,1./Gains)),':');hold on;
Tmpa=get(gca,'YTick');
Tmp=get(gca,'YTickLabel');
set(gca,'YTickLabel',num2str(2.^str2num(Tmp),'%2.2f'));
set(gca,'YTick',Tmpa);
subplot(2,2,3);
plot(T1acc,'.');hold on;
plot(MinRMSI,T1acc(MinRMSI),'r*');
subplot(2,2,4);
plot(T1acc,RMS,'*');hold on;
plot(T1acc(MinRMSI),RMS(MinRMSI),'r*');