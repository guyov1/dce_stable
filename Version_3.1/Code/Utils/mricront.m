function mricront(M,M2,CLRs,Range)
if(nargin<3)
    CLRs=[0 1];
end
if(~exist('M2','var'))
    M2=[];
    if(~exist('Range','var'))
        Range=cell(1,1);
    end
else
    if(~exist('Range','var'))
        Range=cell(1,2);
    end
end
if(~exist([getComputerParams('temppath') filesep 'XX.nii'],'file'))
    copyfile(getKthElement(getComputerParams('tpm'),1),[getComputerParams('temppath') filesep 'XX.nii']);
end
Raw2Nii(M,[getComputerParams('temppath') filesep 'XX.nii'],'float32');
if(~isempty(M2))
    Raw2Nii(M2,[getComputerParams('temppath') filesep 'XX2.nii'],'float32');
    mricronx({[getComputerParams('temppath') filesep 'XX.nii'],[getComputerParams('temppath') filesep 'XX2.nii']},CLRs,Range);
else
    mricronx([getComputerParams('temppath') filesep 'XX.nii'],CLRs,Range);
end