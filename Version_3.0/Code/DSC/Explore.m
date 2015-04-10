function Explore(APath)
if(nargin==0)
    Explore(pwd);
    return;
end
if(isstruct(APath))
    APath=APath.Path;
end
% eval(['!explorer "' APath '" &']);
if(exist(APath,'dir'))
    if(filesep=='/')
    system(['nautilus "' APath '" &']);
    else
         system(['explorer.exe "' APath '" &']);
    end
    return;
end
[BPath BLa Bla2]=fileparts(APath);
if(filesep=='/')
    system(['nautilus "' BPath '" &']);
else
    system(['explorer.exe "' BPath '" &']);
end