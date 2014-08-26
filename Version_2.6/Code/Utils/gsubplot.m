% function gsubplot(N,k)
function gsubplot(N,k,b,c)
if(nargin==2)
    C=ceil(sqrt(N));
    R=ceil(N/C);
    subplot(R,C,k);
else
    subplot(N,k,(b-1)*k+c);
end