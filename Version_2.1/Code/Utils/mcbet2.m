function [TargetFN, MaskTargetFN]=mcbet2(Src,Force,Thresh)
% function [TargetFN MaskTargetFN]=mcbet(Src,Force)
if(~exist('Force','var'))
    Force=0;
end
if(~exist('Thresh','var'))
    Thresh=0.5;
end

[PATHSTR,SName Ext]=fileparts(Src);
TargetFN=[PATHSTR filesep 'BetStripped_' SName Ext];
MaskTargetFN=[PATHSTR filesep 'BetMask_' SName Ext];

if(exist(TargetFN,'file') && ~Force)
    return;
end
delete(TargetFN);
delete(MaskTargetFN);
delete([TargetFN '.gz']);
if(filesep=='/') % Unix
    system(['bet ' Src ' ' TargetFN]);
    system(['gzip -d ' TargetFN '.gz']);
    if(~exist(TargetFN))
        fid=fopen('tmp.sh','w');
        fprintf(fid,'FSLDIR=/usr/share/fsl/5.0\n');
        fprintf(fid,'. ${FSLDIR}/etc/fslconf/fsl.sh\n');
        fprintf(fid,['/usr/share/fsl/5.0/bin/bet ' Src ' ' TargetFN ' -f ' num2str(Thresh) '\n']);
        fclose(fid);
        system('chmod +x tmp.sh');
        system('./tmp.sh');
        delete('tmp.sh');
        system(['gzip -d ' TargetFN '.gz']);
    end
else % Windows
    [A B C]=loadniidata(Src);
    C.untouch=0;
    saveniidata(C,[Src(1:end-4) '.img']);
    system([getComputerParams('mricronpath') 'bet "' Src(1:end-4) '" "' TargetFN(1:end-4) '" -f ' num2str(Thresh)]);
    delete([Src(1:end-4) '.img']);
    delete([Src(1:end-4) '.hdr']);
    [A B C]=loadniidata([TargetFN(1:end-4) '.hdr']);
    Raw2Nii(A,TargetFN,'float32',Src);
    [A B C]=fileparts(Src);
    delete([A filesep 'BetStripped_' B '.img']);
    delete([A filesep 'BetStripped_' B '.hdr']);
end

Mask=loadniidata(TargetFN);
[Vol Hdr A]=loadniidata(Src);
A.img(Mask==0)=0;
saveniidata(A,TargetFN);
A.img(Mask~=0)=1;
saveniidata(A,MaskTargetFN);