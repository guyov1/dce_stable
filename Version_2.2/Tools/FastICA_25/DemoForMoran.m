[sig,mixedsig]=demosig();
figure(1);clf;
for i=1:4
    subplot(4,1,i);
    plot(sig(i,:));
end

figure(2);clf;
for i=1:4
    subplot(4,1,i);
    plot(mixedsig(i,:));
end

[Out1, Out2, Out3] = fastica(mixedsig);

figure(3);clf;
for i=1:4
    subplot(4,1,i);
    plot(Out1(i,:));
end
