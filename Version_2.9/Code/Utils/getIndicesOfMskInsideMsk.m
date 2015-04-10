function Idx=getIndicesOfMskInsideMsk(Small,Big)
Small=Small & Big;
Fs=find(Small);
Fg=find(Big);
Idx=find(ismember(Fg,Fs));