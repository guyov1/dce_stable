function [FA FAPlus FAMinus GoodDiscriminantaPlus GoodDiscriminantaMinus]=CalcFAgSigT1M0(Sig,T1,M0,TR,RoundV,MinV,MaxV)
%               -TR
%                --
%                T1
%     Pd * (1 - e  ) * sin(a)
%     -----------------
%                   -TR
%                    --
%                    T1
%       1  - cos(a)*e

% S=M0(1-E1)SinFA/(1-E1CosFA)
% SinFA is x, CosFA is sqrt(1-x^2), alpha is M0(1-E1)
% S(1-E1(1-tt)/(1+tt))=alpha*2t/(1+tt) : /S
% 1-E1(1-tt)/(1+tt)=beta*t/(1+tt) :*(1+tt), beta is alpha*2/S
% 1+tt-E1+E1tt=beta*t
% tt(1+E1)+t(-beta)+1-E1=0 :
E1=exp(-TR./T1);
a=1+E1;
b=-M0.*(1-E1).*2./Sig;
c=1-E1;
discriminanta=sqrt(b.^2-4.*a.*c);
tPlus=(-b+discriminanta)./(2.*a);
tMinus=(-b-discriminanta)./(2.*a);
GoodDiscriminantaPlus=(imag(tPlus)==0) & isfinite(tPlus);
GoodDiscriminantaMinus=(imag(tMinus)==0)& isfinite(tMinus);
% tPlus=tPlus(imag(tPlus)==0);
% tMinus=tMinus(imag(tMinus)==0);
FAPlus=asind(2.*tPlus(GoodDiscriminantaPlus)./(1+tPlus(GoodDiscriminantaPlus).^2));
FAMinus=asind(2.*tMinus(GoodDiscriminantaMinus)./(1+tMinus(GoodDiscriminantaMinus).^2));
AllFA=[FAPlus; FAMinus];
AllFA=AllFA(imag(AllFA)==0);
Tmp=round(AllFA*RoundV)/RoundV;
FA=mode(Tmp(Tmp>MinV & Tmp<MaxV));