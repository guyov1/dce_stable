
%____________________________________________________________________
%DCE
cd \\fmri-guy2\Dropbox\University\Msc\Thesis\SourceForge\Development
%AddNewScans
%AddNewScans('/mnt/users/Moran/DCE/BVZ_Project/11_Magor_Ilan/W16')
DCEInit
setComputerParamM('temppath','D:\Temp\')
Tmp=getComputerParams('tpm');
copyfile(Tmp{1},[getComputerParams('temppath') 'XX.nii']);
dbstop if error
MainGUI

%rm -rf /mnt/users/Moran/dcescripts/Stable_Versions/code/Version_1.3/dafna
%%%- rm PathDefForDCE folder
%____________________________________________________________________

%____________________________________________________________________
%MRS Single voxel
% cd /data/Moran/MRS_loc
% %Axis controls the 3d order, supposed to be automatic but we can change that in case of mistakes
% %____________________________________________________________________
% 
% MRS Multi voxel
% 
% GUI  /data/Gilad/MRS/MRS_Gui
% Full test data:  /data/Gilad/MRS/SV_CSI
% 
% good localizer: /data/home/gilad/MRS_Loc_For_CSI
% %____________________________________________________________________
