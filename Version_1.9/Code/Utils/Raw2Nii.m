% function Raw2Nii(FN,DIM,Format, [HdrFile])
% Or function Raw2Nii(Data,FNOUT,Format, [HdrFile])
% ---
% the data is saved as "RAW" data. the original DWI are 16 bits integer.
% the processed data (like FA etc.) are 32 bits float. 
function Raw2Nii(FN,DIM,Format, HdrFile)

% FN='D:\Projects\Development\temp\J-0-1.5_AIR.fa';
if(ischar(DIM))
    FNOUT=DIM;
    I=FN;
    DIM=size(I);
else
    [A B C]=fileparts(FN);
    FNOUT=fullfile(A,[B '.nii']);
    F=fopen(FN,'r');
    [Data N]=fread(F,Inf,'float');
    I=reshape(Data, DIM);
    I=flipdim(flipdim(I,2),1);
end
if(nargin<4)
    Base=getComputerParams('temppath');
%     HdrFile='/u/peptibase3-ext/libermg1/Temp/ForT.nii';
%     HdrFile=[Base '/ForT.nii'];
%     HdrFile=['/data/Gilad/ForT.nii'];
%     HdrFile=['C:\STRDCE\John\Database\DCEOut\KaRo_20070911\Baseline.nii'];
    HdrFile=[getComputerParams('temppath') filesep 'XX.nii'];
end
[A B NII]=loadniidata(HdrFile);
NII.hdr.dime.dim=[numel(DIM) DIM ones(1,7-numel(DIM))];
NII.hdr.dime.scl_slope=1;
switch(Format) % see help save_nii
    case 'float32'
        NII.hdr.dime.datatype=16;
        NII.hdr.dime.bitpix=32;
    case 'int16'
        NII.hdr.dime.datatype=4;
        NII.hdr.dime.bitpix=16;
    case 'uint16'
        NII.hdr.dime.datatype=512;
        NII.hdr.dime.bitpix=16;
    case 'uint8'
        NII.hdr.dime.datatype=2;
        NII.hdr.dime.bitpix=8;
    otherwise
        error('Unknown bit precision');
end
NII.img=I;
saveniidata(NII,FNOUT);
