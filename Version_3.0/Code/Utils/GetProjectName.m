function Out=GetProjectName(Path)
% S=regexp(Path,'\','split');
S=regexp(Path,filesep);
N=getComputerParams('ProjectNamePathDepth');
if(length(S)>N)
    Out=Path((S(N)+1):(S(N+1)-1));
else
    Out='None';
end
% Out=S{3};