function [MnV MV]=FindDR(Data)
FData=Data(isfinite(Data));
FData=FData(FData~=0);
FData=FData(FData~=min(FData));
FData=FData(FData~=max(FData));
FData=getKrandomSamples(FData,min(numel(FData),10000));
DR=max(FData)-min(FData);
NData=(FData-min(FData))/DR;

UVals=unique(NData);
NUVals=(UVals-min(UVals))./(max(UVals)-min(UVals));
for i=1:numel(UVals); CDF(i)=sum(NData(:)<=UVals(i));end
NCDF=CDF./CDF(end);
W=1;
Energy=W*NCDF-NUVals';
[F FIx]=max(Energy(2:end));
MV=UVals(FIx+1)*DR+min(FData);

[F FIn]=min(Energy(1:FIx-1));
MnV=UVals(FIn)*DR+min(FData);
