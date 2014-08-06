function CTC=DCEtoCTCf(Msk,DCE4D,CT1,FA3D,DCETR, Baseline)

R10=1./CT1;
E1=exp(-DCETR.*R10);

CosFA=cosd(FA3D);
SinFA=sind(FA3D);

DCEM0=(Baseline.*(1-E1.*CosFA))./((1-E1).*SinFA);
M0Clmn=DCEM0(Msk);

DCE2D=Reshape4d22d(DCE4D,Msk);

TmpVar=repMulti(DCE2D,1./M0Clmn);
E1F=(repPlus(TmpVar,-SinFA(Msk)))./(repPlus(repMulti(TmpVar,CosFA(Msk)),-SinFA(Msk)));
R1F=-log(E1F)./DCETR;
R10Clmn=R10(Msk);
CTC=repPlus(R1F,-R10Clmn);