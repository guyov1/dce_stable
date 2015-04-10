disp('Updating very short infos');
Q=load(getComputerParams('ShortInfosFN'));
Arr=structFNs2structArr(Q.ShortInfos);
N={Arr.Name};
Names=unique(N);
save(getComputerParams('VeryShortInfosFN'),'Names');