function Out=gDicom2Nifti(SRC_DIR,TRG_FN)
% Return values:
% 1 - Already exists
% 

TmpName= license('inuse');
TmpName=TmpName(1).user;
Base=[getComputerParams('temppath') '_DCM2Nii_' TmpName '_' regexprep(datestr(now),'[-: ]','_') '_' num2str(floor(rand*1000)) filesep];
mkdir(Base);

if (filesep ~='/') % If not Unix
    if(TRG_FN(1)==filesep  && TRG_FN(2) ~=filesep)
        TRG_FN=[filesep TRG_FN];
    end
    if(SRC_DIR(1)==filesep  && SRC_DIR(2) ~=filesep)
        SRC_DIR=[filesep SRC_DIR];
    end
end

if(exist(TRG_FN,'file'))
    disp(['Already exist ' TRG_FN]);
    Out=1;
else
    try
        jkhkjh;
        Out=mricroDcm2nii(SRC_DIR,TRG_FN);
        Out=numel(dir([TRG_FN(1:end-4) '*']));
        %        disp(['No nifti for' SRC_DIR]);
    catch Error_Struct
        s = warning('query', 'all');
        delete(fullfile(Base,'func/run_0001/',['*.nii']));
        warning('off','MATLAB:FINITE:obsoleteFunction');
        Out=dicom2niftix('dicom_dir',SRC_DIR,'subject_dir',Base);
        warning(s);

        if(Out>0)
            if(exist(fullfile(Base,'anat','overlay.nii'),'file'))
                movefile(fullfile(Base,'anat','overlay.nii'),TRG_FN);
            else
                if(exist(fullfile(Base,'anat','anatomy.nii'),'file'))
                    movefile(fullfile(Base,'anat','anatomy.nii'),TRG_FN);
                else
                    if(exist(fullfile(Base,'dti','dti_01.nii'),'file') || exist(fullfile(Base,'dti','dti_0001.nii'),'file'))
                        D=dir(fullfile(Base,'dti','dti_*.nii'));
                        NFiles=length(D);
                        for i=1:NFiles
                            Suffix=D(i).name(4:end);
                            [TPath TName TExt]=fileparts(TRG_FN);
                            CurTFN=fullfile(TPath,[TName Suffix]);
                            movefile([Base filesep 'dti' filesep D(i).name],CurTFN);
                        end
                    else
                        [A B C]=fileparts(TRG_FN);
                        D=dir(fullfile(Base,'func/run_0001/',[B '_*.nii']));
                        if(isempty(D))
                            D=dir(fullfile(Base,'func/run_0001/',['vol' '_*.nii']));
                        end
                        if(~isempty(D))
                            NFiles=length(D);
                            for i=1:NFiles
                                Suffix=D(i).name(4:end);
                                [TPath TName TExt]=fileparts(TRG_FN);
                                CurTFN=fullfile(TPath,[TName Suffix]);
                                movefile([Base filesep 'func/run_0001' filesep D(i).name],CurTFN);
                            end
                        else
                            D=dir(Base);
                            D=D(3:end);
                            D=D([D.datenum]>(now-5/(24*60)));
                            D=D(~ismember({D.name},{'Coreg' 'anat' 'dti' 'func' 'N3'}));
                            if(numel(D)==1)
                                DD=dir([Base filesep D.name filesep]);
                                movefile([Base filesep D.name filesep DD(3).name],TRG_FN);
                                delete([Base filesep D.name filesep '*.*']);
                                rmdir([Base filesep D.name]);
                            else
                                disp(['Error no file created - ' SRC_DIR ' - ' TRG_FN]);
                                error(['Error no file created - ' SRC_DIR ' - ' TRG_FN]);
                            end
                        end
                    end
                end
            end
        end
    end
end
pause(0.5);
rmdir(Base,'s');