function FoundFN=getNiiFNFromInfo(Info)
N=numel(Info);
if(N>1)
    FoundFN=cell(1,N);
    for i=1:N
        FoundFN{i}=getNiiFNFromInfo(getKthElement(Info,i));
    end
    FoundFN=flattenCell(FoundFN);
else
    if(~strcmp(Info.Class,'SWI'))
        FoundFN=getNiiFullFNFromFN(Info.SeriesInstanceUID);
    else
        FoundFN=getNiiFullFNFromFN([Info.SeriesInstanceUID '_' sprintf('%02d',4*(Info.EchoNumber-1)+1)]);
%         FoundFN=getNiiFullFNFromFN(Info.SeriesInstanceUID,1);
        if(isempty(FoundFN)) % couldn't find mutli nifti SWI
            FoundFN=getNiiFullFNFromFN(Info.SeriesInstanceUID);
            if(~isempty(FoundFN))
                disp(['getNiiFNFromInfo - SWI without phase! ' Info.Path]);
            end
        end
    end
end