function [Cost ParamVec]=FindB1gWMCoeffsWithNonLinTransform(UncleanedData,ProbMat,W,NonLinTransform,rNonLinTransform)
SolVec=NonLinTransform(UncleanedData);
ParamVec=(repMulti(ProbMat,W))\(W.*SolVec);
Errs=rNonLinTransform(SolVec)-rNonLinTransform(ProbMat*ParamVec);
Cost=sqrt((Errs'.^2)*(W/sum(W)));