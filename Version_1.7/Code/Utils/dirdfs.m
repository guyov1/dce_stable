function [Out N]=dirdfs(BaseDir,Prefix,h,N)
Base=0;
if(nargin==1)
    Prefix='';
    N=0;
    h = waitbar(0,num2str(0));
    Base=1;
end
DD=dir(BaseDir);
DD=DD(3:end);
Out=struct([]);
waitbar(0.5,h,num2str(N)); 
for i=1:length(DD)
    NameNoPrefix=DD(i).name;
    DD(i).name=[Prefix DD(i).name];
    Out=[Out; DD(i)];
    N=N+1;
    waitbar(0.5,h,num2str(N));
    if(DD(i).isdir)
        [CurAdd N]=dirdfs([BaseDir filesep NameNoPrefix], [DD(i).name filesep],h,N);
%         for j=1:length(CurAdd)
            Out=[Out; CurAdd];
%         end
    end
end
if(Base)
    close(h);
end
