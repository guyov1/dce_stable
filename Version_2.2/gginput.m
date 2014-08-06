function Out=gginput
but=1;
Out=zeros(0,2);
while(but==1)
    [x,y,but] = ginput(1);
    if(numel(x)>1)
        but=3;
    end
    Out=[Out; [x y]];
end