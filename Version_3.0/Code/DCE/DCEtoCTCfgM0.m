function CTC=DCEtoCTCfgM0(Msk,DCE4D,CT1,FA3D,DCETR, DCEM0)

R10=1./CT1;
% E1=exp(-DCETR.*R10);

CosFA=cosd(FA3D);
SinFA=sind(FA3D);

% Get the masked M0 values
% DCEM0=(Baseline.*(1-E1.*CosFA))./((1-E1).*SinFA);
M0Clmn=DCEM0(Msk);

% Reshape the 4d image to 2d
DCE2D=Reshape4d22d(DCE4D,Msk);

%S(t)/ M0 = sin(FA)*(1-E1(t))/(1-E1(t)cos(FA))
TmpVar=repMulti(DCE2D,1./M0Clmn);

% ( 1-E1(t)cos(FA) ) *  S(t)/ M0 = sin(FA)*(1-E1(t))    =>
%  E1(t) * ( sin(FA) -  cos(FA) * S(t)/ M0 ) = sin(FA) - S(t)/ M0 )  =>
%  E1(t) = ( sin(FA) - S(t)/ M0 ) /  (  sin(FA) -  cos(FA) * S(t)/ M0 ) =>
%  E1(t) = ( -sin(FA) + S(t)/ M0 ) /  ( -sin(FA) +  cos(FA) * S(t)/ M0 )
E1F=(repPlus(TmpVar,-SinFA(Msk)))./(repPlus(repMulti(TmpVar,CosFA(Msk)),-SinFA(Msk)));

% E1(t) = exp ( -TR/T1(t) ) ->   1 / T1(t) = -log(E1) / TR  , R1(t) = 1 / T1(t)
R1F=-log(E1F)./DCETR;

% R0 = 1 / T0
R10Clmn=R10(Msk);

% 1/T1(t)=1/(T1base)+C(t) =>
% C(t) = 1/T1(t) - 1/(T1base)
CTC=repPlus(R1F,-R10Clmn);