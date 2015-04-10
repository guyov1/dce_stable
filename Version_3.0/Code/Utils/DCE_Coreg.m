DCECoregP=[WorkingP 'DCEMainCoreged' filesep];
mkdir(DCECoregP);
DDCE=dir([DCEMNiiOrigP '*.nii']);
DCEFNs=strcat(DCEMNiiOrigP,{DDCE.name})';
DCEFNs=DCEFNs((Options.nVolsToRemoveFromStart+1):end);
if(Options.MainCoregistration==1)
        MatFN=RealignEstimate(DCEFNs,Force,false);
        DCECrgFNs=CoregWrite(DCEFNs,MatFN,Force,DCECoregP,false);
end
if(Options.MainCoregistration==0)
        OrigP=[WorkingP 'DCEMainNii' filesep];
        TrgP=[WorkingP 'DCEMainCoreged' filesep];
        D=dir([OrigP '*.nii']);
        for i=1:numel(D)
            copyfile([OrigP D(i).name],[TrgP 'Coreged_' D(i).name]);
        end
end
if(Options.MainCoregistration>1)
    OrigP=[WorkingP 'DCEMainNii' filesep];
    TrgP=[WorkingP 'DCEMainCoreged'];
    D=dir([OrigP '*.nii']);
    TrgVolFN=[OrigP D(Options.MainCoregistration).name];
    for i=1:numel(D)
        disp(['Coregstering ' num2str(i) ' out of ' num2str(numel(D))]);
        MainMatFNs{i}=CoregEstimate([OrigP D(i).name],TrgVolFN,Force);
        MainCrgFNs{i}=CoregWrite([OrigP D(i).name],MainMatFNs{i},Force,TrgP,false,TrgVolFN);
    end
end