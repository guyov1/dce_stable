% This functions is being called from ReadNewScans.m or UpdateShortInfos.m
function Out=Infos2ShortInfos(Infos)

FNs=fieldnames(Infos);
ShortFields={'SeriesInstanceUID','Filename','FileModDate','StudyDate','SeriesDate','StudyTime','SeriesTime','SeriesDescription','PatientName','PatientAge'};
ShortFields=[ShortFields 'RepetitionTime','EchoTime','InversionTime','EchoNumber','EchoTrainLength'];
ShortFields=[ShortFields 'FlipAngle','SeriesNumber','ImagesInAcquisition'];

for i=1:length(FNs)
    for j=1:length(ShortFields)
        if(isfield(Infos.(FNs{i}),ShortFields{j}))
            CurShort.(ShortFields{j})=Infos.(FNs{i}).(ShortFields{j});
        else
            CurShort.(ShortFields{j})=-777;
        end
    end
    if(isfield(Infos.(FNs{i}),'Private_0019_10f9'))
        CurShort.TG=str2num(char((Infos.(FNs{i}).Private_0019_10f9)'));
    else
        CurShort.TG=-777;
    end
    % RAS
    CurShort.Position='None';
    if(isfield(Infos.(FNs{i}),'Private_0019_1018'))
        switch(Infos.(FNs{i}).Private_0019_1018(1))
            case {73,83} % 'I,S - Ax
                CurShort.Position='AX';
             case {76,82} % R,L - Sag
                CurShort.Position='SAG';
            case {65,80} % A,P - Cor
                CurShort.Position='COR';
            otherwise
                CurShort.Position='None';
        end
    end
    
    CurShort.ImagesInAcquisition=-777;
    if(isfield(Infos.(FNs{i}),'Private_0025_1007'))
%         CurShort.ImagesInAcquisition=Infos.(FNs{i}).Private_0025_1007(1);
        if(numel(Infos.(FNs{i}).Private_0025_1007')==1)
            CurShort.ImagesInAcquisition=Infos.(FNs{i}).Private_0025_1007;
        else
            CurShort.ImagesInAcquisition=double(Infos.(FNs{i}).Private_0025_1007')*(256.^[0:3]');
        end
        if(numel(CurShort.ImagesInAcquisition)>1)
            aq=232;
        end
    end
    if(isfield(Infos.(FNs{i}),'Private_2001_1018'))
        CurShort.ImagesInAcquisition=double(Infos.(FNs{i}).Private_2001_1018);
    else
        a=5;
    end
    CurShort.R1=-777;
    if(isfield(Infos.(FNs{i}),'Private_0019_108a'))
        CurShort.R1=Infos.(FNs{i}).Private_0019_108a(1);
    end
    CurShort.R2=-777;
    if(isfield(Infos.(FNs{i}),'Private_0019_108b'))
        CurShort.R2=Infos.(FNs{i}).Private_0019_108b(1);
    end
    Out.(FNs{i})=CurShort;
end