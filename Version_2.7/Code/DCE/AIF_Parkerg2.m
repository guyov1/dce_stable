function C=AIF_Parkerg2(t,A1,sig1,T1,rA2,sig2,T2,ralpha,beta)

% 	A_1    � First Gaussian amplitude
% 	T_1    � First Gaussian delay
% 	sig_1  � First Gaussian lobe width 
% 	A_2    � Second Gaussian amplitude
% 	T_2    � Second Gaussian delay
% 	sig_2  � Second Gaussian lobe width 
% 	alpha  � exponent amplitude (starting decay value)
% 	T      - sigmoid shift (starting point for the decay procedure)
% 	beta   � exponent decay rate
% 	s      � sigmoid decay rate

% Latex equation
%\[{C_p}\left( t \right) = \mathop \sum \limits_{n = 1}^2 \frac{{{A_n}}}{{{\sigma _n}{{\left( {2\pi } \right)}^{\frac{1}{2}}}}}{e^{\frac{{ - {{\left( {t - {T_n}} \right)}^2}}}{{2{\sigma _n}^2}}}} + \alpha \frac{{{e^{ - \beta t}}}}{{1 + {e^{ - s\left( {t - \tau } \right)}}}}\]


A2    = rA2*A1;
alpha = ralpha*A1;

% http://dx.doi.org/10.1002/mrm.21066
% A1=0.809;A2=0.330;T1=0.17046;T2=0.365;
% sig1=0.0563;sig2=0.132;alpha=1.050;beta=0.1685;
% s=38.078;tau=0.483;
% A~mM.min, T~min; sig~min, alpha~mM, beta,s~1/min, tau~min

s2pi = sqrt(2*pi);
C1   = A1.*(sig1*s2pi).*exp(-((t-T1).^2)./(2*(sig1^2)))./(sig1*s2pi);
C2   = A2.*(sig2*s2pi).*exp(-((t-T2).^2)./(2*(sig2^2)))./(sig2*s2pi);

% MaxAmp=max(C);
% MaxAmp=A1/(sig1*s2pi);

% SigmoidPart1=exp(-beta*t)./(1+exp(-s*(t-tau)));
SigmoidPart = (t>T1).*exp(-beta*(t-T1));
% MaxSigmoidPart=max(SigmoidPart);

% By taking the max, Gilad wanted to avoid changing
% of one gaussian paramters to affect the second one
C = max(C1,C2+alpha*SigmoidPart);

% MaxSigmoidPart=1;
% C=C + alpha*SigmoidPart;

% C=max(A1.*exp(-((t-T1).^2)./(2*(sig1^2)))./(sig1*s2pi),A2.*exp(-((t-T2).^2)./(2*(sig2^2)))./(sig2*s2pi))+alpha*exp(-beta*t)./(1+exp(-s*(t-tau)));

% C=max(max(A1.*exp(-((t-T1).^2)./(2*(sig1^2)))./(sig1*s2pi),A2.*exp(-((t-T
% 2).^2)./(2*(sig2^2)))./(sig2*s2pi)),alpha*exp(-beta*t)./(1+exp(-s*(t-tau)
% )));