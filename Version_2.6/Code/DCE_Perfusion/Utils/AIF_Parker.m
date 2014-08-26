function Ca = AIF_Parker(t,A1,sig1,T1,A2,sig2,T2,alpha,beta,s,tau,additional_time_shift)

% http://dx.doi.org/10.1002/mrm.21066
% A1=0.809;A2=0.330;T1=0.17046;T2=0.365;
% sig1=0.0563;sig2=0.132;alpha=1.050;beta=0.1685;
% s=38.078;tau=0.483;

% A~mM.min, T~min; sig~min, alpha~mM, beta,s~1/min, tau~min

s2pi = sqrt(2*pi);

num_AIFs        = length(additional_time_shift);
num_time_points = length(t);

% Check if duplication is needed
if  num_AIFs> 1
    
    % Create relevant matrices
    dup_t     = repmat(t,num_AIFs,1);
    dup_shift = transpose(repmat(additional_time_shift,num_time_points,1));
    
    Ca    = A1.*exp(-(((dup_t-dup_shift)-T1).^2)./(2*(sig1^2)))./(sig1*s2pi) + ...
            A2.*exp(-(((dup_t-dup_shift)-T2).^2)./(2*(sig2^2)))./(sig2*s2pi) + ...
            alpha*exp(-beta*(dup_t-dup_shift))./(1+exp(-s*((dup_t-dup_shift)-tau)));
else
    Ca    = A1.*exp(-(((t-additional_time_shift)-T1).^2)./(2*(sig1^2)))./(sig1*s2pi) + ...
            A2.*exp(-(((t-additional_time_shift)-T2).^2)./(2*(sig2^2)))./(sig2*s2pi) + ...
            alpha*exp(-beta*(t-additional_time_shift))./(1+exp(-s*((t-additional_time_shift)-tau)));
end



% C=max(A1.*exp(-((t-T1).^2)./(2*(sig1^2)))./(sig1*s2pi),A2.*exp(-((t-T2).^2)./(2*(sig2^2)))./(sig2*s2pi))+alpha*exp(-beta*t)./(1+exp(-s*(t-tau)));

% C=max(max(A1.*exp(-((t-T1).^2)./(2*(sig1^2)))./(sig1*s2pi),A2.*exp(-((t-T2).^2)./(2*(sig2^2)))./(sig2*s2pi)),alpha*exp(-beta*t)./(1+exp(-s*(t-tau))));