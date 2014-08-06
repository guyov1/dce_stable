function ApproximatePeakTime=FindApproximatePeakTime(CTC2D,BolusStart,NumVols)

[MaxVal, PeakTime]=max(CTC2D,[],2);
GoodPeak=PeakTime>(BolusStart-5) & PeakTime<(BolusStart+15);
BeforeVals=[CTC2D(sub2ind(size(CTC2D),find(GoodPeak),max(1,PeakTime(GoodPeak)-1)))];
AfterVals=[CTC2D(sub2ind(size(CTC2D),find(GoodPeak),min(NumVols,PeakTime(GoodPeak)+1)))];
NearVals=[BeforeVals AfterVals];
MaxVal1D=MaxVal(GoodPeak);
PeakTime1D=PeakTime(GoodPeak);
[MaxNearVal, WhichDirection]=max(NearVals,[],2);
RelativeDist=MaxNearVal./(MaxVal1D+MaxNearVal);
ApproximatePeakTime1D=PeakTime(GoodPeak)+((WhichDirection-1)*2-1).*RelativeDist;
ApproximatePeakTime=PeakTime;
ApproximatePeakTime(GoodPeak)=ApproximatePeakTime1D;
ApproximatePeakTime(~GoodPeak)=NaN;
