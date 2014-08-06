%function mricronx( in , Color, Range)
function mricronx( Cin , Color, Range)
if(~iscell(Cin))
    Cin={Cin};
end
N=numel(Cin);
if(~exist('Color','var'))
    Color=0;
end
if(numel(Color)==1)
    Color=repmat(Color,[N 1]);
end
if(~exist('Range','var'))
    Range=cell(1,N);
end
if(~iscell(Range))
    Range={Range};
end
if(isempty(Range)==1)
    Range=cell(1,N);
end
FNs=cell(1,N);
CLR=cell(1,N);
RNG=cell(1,N);
for i=1:numel(Cin)
    in=Cin{i};
    if(isstruct(in))
        FNs{i}=getNiiFNFromInfo(in);
        Z=true;
    else
        FNs{i}=in;
        Z=~isempty(dir(in));
    end
    if(isempty(FNs{i}) || ~Z)
        disp('mricron - file not found');
        return;
    end
    CLR{i}=[' -c -' num2str(Color(i))];
    if(~isempty(Range{i}))
        RNG{i}=[' -l ' num2str(Range{i}(1)) ' -h ' num2str(Range{i}(2))];
    else
        RNG{i}=' ';
    end
end

if(filesep=='/') % Unix
     CallStr=[getComputerParams('mricronpath') 'mricron'];
     CallStr=['mricron'];
else % windows
    CallStr=['"' getComputerParams('mricronpath') 'mricron.exe"'];
end

A=[CallStr ' "' FNs{1} '" ' CLR{1} RNG{1}];
% -b is % 0 20 40 60 80 100 or -1 additive blending
for i=2:N
    A=[A ' -o "' FNs{i} '" ' CLR{i} RNG{i} ' -b -1'];
end
A=[A ' &'];
% disp(A);
system(A);