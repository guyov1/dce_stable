function Out=ClassifySeq(Path,X)
Out='';
if(length(Path)<2)
    return;
end
[base, last]=fileparts(Path);
%      if(strhas(last,'Series'))
%          disp(['Needs renaming - ' last]);
%      end
if(nargin<2 && (~(strhas(last,'Se') && ~(strhas(last,'_Case') || strhas(last,'Baseline')) ) ) )
    Out='Klum';
else
    if(strhas(last,'FIESTA'))
        if(~isempty(Out))
            disp(['Error - ' last]);
        end
        Out='FIESTA';
    end
    if(strhas(last,'MRS'))
        if(~isempty(Out))
            disp(['Error - ' last]);
        end
        Out='MRS';
    end
    if(strhas(last,' O2'))
        if(~isempty(Out))
            disp(['Error - ' last]);
        end
        Out='O2';
    end
    if(strhas(last,'CO2'))
        if(~isempty(Out))
            disp(['Error - ' last]);
        end
        Out='CO2';
    end
    if(strhas(last,'DCE') || strhas(last,'T1 with FA') || strhas(last,'Uptake') || strhas(last,'DESPOT'))
        if(~isempty(Out))
            disp(['Error - ' last]);
        end
        Out='DCE';
    end
    if(strhas(last,'T1MAP'))
        if(~isempty(Out))
            disp(['Error - ' last]);
        end
        Out='T1MAP';
    end
    if((strhas(last,'SWI') && ~strhas(last,'T2 FSE')) || strhas(last,'cctrl') || strhas(last,'cclrl'))
        if(~isempty(Out))
            disp(['Error' last]);
        end
        Out='SWI';
    end
    if(strhas(last,'DSC') || strhas(last,'DSE'))
        if(~isempty(Out))
            disp(['Error' last]);
        end
        Out='DSC';
    end
    if(strhas(last,'FLAIR'))
        if(~isempty(Out))
            disp(['Error' last]);
        end
        Out='FLAIR';
    end
    if(strhas(last,'spgr') || strhas(last,'Nav'))
        if(~isempty(Out))
            disp(['Error' last]);
        end
        Out='T13DSPGR';
    end
    if((strhas(last,'T1 ') || strhas(last,' T1') || strhas(last,'flip')) && ~strhas(last,'DCE') && ~strhas(last,'spgr') && ~strhas(last,'T1 with FA') && ~strhas(last,'DESPOT1') && ~strhas(last,'Uptake'))
        if(~isempty(Out))
            disp(['Error' last]);
        end
        Out='T1';
    end
    if(strhas(last,'LOC'))
        if(~isempty(Out))
            disp(['Error' last]);
        end
        Out='LOC';
    end
    if(strhas(last,'DTI'))
        if(~isempty(Out))
            disp(['Error' last]);
        end
        Out='DTI';
    end
    if(strhas(last,'T2 FSE'))
        if(~isempty(Out))
            disp(['Error' last]);
        end
        Out='T2_FSE';
    end
    if(strhas(last,'T2 \+PD'))
        if(~isempty(Out))
            disp(['Error' last]);
        end
        Out='T2_Plus_PD';
    end
    if(strhas(last,'SCREENSAVE'))
        if(~isempty(Out))
            disp(['Error' last]);
        end
        Out='SCREENSAVE';
    end
    if(strhas(last,'fMRI-RT'))
        if(~isempty(Out))
            disp(['Error' last]);
        end
        Out='fMRI_RT';
    end
    if(strhas(base,'Analysis'))
        if(~isempty(Out) && ~(strhas(last,'DCE') || strhas(last,'DSC-Perfusion')))
            disp(['Error' txt{i}]);
        end
        Out='Analysis';
    end
end