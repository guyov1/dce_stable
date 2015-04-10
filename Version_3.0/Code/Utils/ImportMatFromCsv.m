FN='Export.csv';

[A, delimiter]=importdata(FN);
B=importdata(FN,delimiter,size(A.data,1)+2);
Export.CurPKs=A.data;
Export.SampleTs=B.data(1,:);
Export.CurCTCs=B.data(2:end,:);