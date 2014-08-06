DCECoregP=[WorkingP 'DCEMainCoreged/'];
mkdir(DCECoregP);
DDCE=dir([DCEMNiiOrigP '*.nii']);
DCEFNs=strcat(DCEMNiiOrigP,{DDCE.name})';
MatFN=RealignEstimate(DCEFNs,Force,false);
DCECrgFNs=CoregWrite(DCEFNs,MatFN,Force,DCECoregP,false);