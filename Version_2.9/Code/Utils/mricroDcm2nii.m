function Out=mricroDcm2nii(SRC_DIR,TRG_FN)
% function Out=gDicom2Nifti(SRC_DIR,TRG_FN)
TmpName= license('inuse');
TmpName=TmpName(1).user;
Base=[getComputerParams('temppath') filesep 'Temp_mricroDCM2Nii_' TmpName '_' regexprep(datestr(now),'[-: ]','_') '_' num2str(floor(rand*1000)) filesep];
% Base=getComputerParams('temppath');
mkdir(Base);
% delete('/local_data/Temp/*.*');
if(filesep=='/') % Unix
    CallStr=[getComputerParams('mricronpath') 'dcm2nii -a n -g n -o ' Base ' "' SRC_DIR '"' ' '];
    %system(CallStr,' -echo');
    system(CallStr);
else % windows
    CallStr=['"' getComputerParams('mricronpath') 'dcm2nii.exe"  -a n -g n -o ' '"' Base '"'  ' "' SRC_DIR '"' ' '];
    system(CallStr);
end
% if(any(strcmp(getComputerParams('name'),{'T3','FranceLabG'})))
% %     system(['/local_data/OutsideCode/mricron/dcm2nii -a n -g n -o ' Base ' "' SRC_DIR '"' ' &'],'-echo');
%     system(['/data/Gilad/mricron/dcm2nii -a n -g n -o ' Base ' "' SRC_DIR '"' ' '],'-echo');
% else
%     system(['"C:\Program Files\MRIcroN\dcm2nii.exe" -a n -g n -o ' Base ' "' SRC_DIR '"' '&']);
% end
% pause(1.5);
D=dir([Base filesep '*.nii']);
nD=numel(D);
if(nD==0)
    disp(['Error no file created - ' SRC_DIR ' - ' TRG_FN]);
    Out=[];
    return;
end
Out=nD;
if(nD==1)
    A=loadniidata([Base filesep D(1).name]);
    nVols=size(A,4);
    if(nVols==1)
        movefile([Base filesep D(1).name],TRG_FN);
    else
        try
            AA=load_untouch_nii([Base filesep D(1).name]);
        catch
            AA=load_nii([Base filesep D(1).name]);
        end
        B=AA;
        B.hdr.dime.pixdim(5)=1;
        B.hdr.dime.dim(5)=1;
        [PFN BFN EFN]=fileparts(TRG_FN);
        for v=1:nVols
            B.img=squeeze(AA.img(:,:,:,v));
            CurFN=fullfile(PFN,[BFN '_' num2str(v,'%04d') EFN]);
            %             Raw2Nii(
            saveniidata(B,CurFN);
        end
    end
else
    for i=1:nD
        Suffix=num2str(i,'%02.0f');
        [TPath TName TExt]=fileparts(TRG_FN);
        CurTFN=fullfile(TPath,[TName '_' Suffix TExt]);
        movefile([Base filesep D(i).name],CurTFN);
    end
end
rmdir(Base,'s');