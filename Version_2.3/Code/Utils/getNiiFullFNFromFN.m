function FoundFN=getNiiFullFNFromFN(inFN,DTI)
if(~exist('DTI','var'))
    DTI=0;
end
if(DTI>0) % Multi volume sequence
    [PATHSTR,NAME,EXT]=fileparts(inFN);
    if(strcmp(EXT,'.nii'))
        FirstFN=[NAME '_01'];
    else
        FirstFN=[NAME EXT '_01'];
    end
    BaseFN=getNiiFullFNFromFN(FirstFN,0);
    if(isempty(BaseFN))
        FoundFN=BaseFN;
        return;
    end
    D=dir([BaseFN(1:end-6) '*.nii']);
    [DPATHSTR]=fileparts(BaseFN);
    FoundFN=cell(numel(D),1);
    for i=1:length(D)
        FoundFN{i}=fullfile(DPATHSTR,D(i).name);
    end
    nSWI=numel(FoundFN);
    if(DTI==1 && nSWI>1) % Only magnitude
        FoundFN=FoundFN(((1:(nSWI/4))*4)-3);
    end
    return;
end
if(exist(which(inFN),'file'))
    FoundFN=which(inFN);
    return;
end
FoundFN='';
[PATHSTR,NAME,EXT]=fileparts(inFN);
for i=1:getComputerParams('nNiftisPaths')
    if(strcmp(EXT,'.nii'))
        FN=fullfile( getComputerParams('NiftisPath',i), [NAME '.nii'] );
    else
        FN=fullfile( getComputerParams('NiftisPath',i), [NAME EXT '.nii'] );
    end
    if(exist(FN,'file'))
        FoundFN=FN;
        break;
    end
end