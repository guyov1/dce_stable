function [ScoreOut FAsOut GainsOut SimOut]=T1M0From2(SOrdered, FAOrdered, TROrdered)
[RelaxBase{1:3}]=CalcT1byFAfw2(SOrdered(:,1:2)',FAOrdered(1:2),TROrdered(1:2));

CurFAs=FAOrdered;
% for Iter=1:1
%     CurSim=SPGRfM(RelaxBase{1}',RelaxBase{2}',CurFAs(3:end),TROrdered(3:end));
%     CurGains=[0 0 median(SOrdered(:,3:end)./CurSim,1)];
%     for ii=3:numel(TROrdered)
%         CurFAs(ii)=CalcFAgSigT1M0(SOrdered(:,ii)/CurGains(ii),RelaxBase{1},RelaxBase{2},TROrdered(ii),100,max(0.5,FAOrdered(ii)/2),FAOrdered(ii)*2);
%     end
% end
for ii=3:numel(TROrdered)
    FARange=max(0.5,FAOrdered(ii)-10):0.1:(FAOrdered(ii)+10);
    FAsSim=SPGRfM(RelaxBase{1}',RelaxBase{2}',FARange,FARange*0+TROrdered(ii));
    clear GainForii
    for j=1:numel(FARange)
        GainForii(j)=FAsSim(:,j)\SOrdered(:,ii);
    end
    [MinSRMS,MIdx]=min(mean(repPlus(repMulti(FAsSim',GainForii'),-SOrdered(:,ii)').^2,2));
    CurFAs(ii)=FARange(MIdx);
    CurGains(ii)=GainForii(MIdx);
end
for ii=3:numel(TROrdered)
    SimOut(:,ii-2)=CurGains(ii)*SPGRfM(RelaxBase{1}',RelaxBase{2}',CurFAs(ii),TROrdered(ii));
end
ScoreOut=sqrt(mean(median((SOrdered(:,3:end)-SimOut).^2,1),2));