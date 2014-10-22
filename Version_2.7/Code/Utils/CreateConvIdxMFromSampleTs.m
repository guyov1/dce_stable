function Out=CreateConvIdxMFromSampleTs(nTimePoints)

[X Y]=meshgrid(1:nTimePoints);
Out=((Y-X)+1).*(X<=Y);
