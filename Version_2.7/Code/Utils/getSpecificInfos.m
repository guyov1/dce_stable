function [Out OutList]=getSpecificInfos(varargin)
% ShortInfos=load(getComputerParams('ShortInfosFN'),'ShortInfos');
% ShortInfos=ShortInfos.ShortInfos;
if(mod(nargin,2)==1)
    Lnargin=nargin-1;
    Lvarargin=varargin(2:end);
    GotInfos=true;
    ShortInfosArr=varargin{1};
else
    GotInfos=false;
    Lnargin=nargin;
    Lvarargin=varargin;
    

    load(getComputerParams('ShortInfosFN'),'ShortInfos');
    
    FNs=fieldnames(ShortInfos);
    ShortInfosArr=ShortInfos.(FNs{1});
    ShortInfosArr(length(FNs))=ShortInfos.(FNs{1});
    for i=1:length(FNs)
        ShortInfosArr(i)=ShortInfos.(FNs{i});
    end
end
if(Lnargin==0)
    Out=ShortInfosArr;
    return;
end
% Cats=fieldnames(ShortInfos.(FNs{1}));
Cats=fieldnames(ShortInfosArr(1));
LCats=lower(Cats);
XCats=MultiConcat('x',LCats);
if(~GotInfos || ~isfield(ShortInfosArr,'Usable'))
    A=ones(1,length(ShortInfosArr))>0;
else
    A=ones(1,length(ShortInfosArr))>0;
    A(~gIsEmpty({ShortInfosArr.Usable}))=[ShortInfosArr.Usable];
%     A=[ShortInfosArr.Usable];
end
for i=1:floor(Lnargin/2)
    CurCat=Lvarargin{i*2-1};
    Arg=Lvarargin{i*2};
    switch lower(CurCat)
        %         case {'Class','StudyDate','SeriesDate','Name','Project'}
        case XCats
            CurCat=CurCat(2:end);
            FieldName=Cats(strcmpi(CurCat,Cats));
            Classes={ShortInfosArr.(FieldName{1})};
            %             B=strcmp(Classes,varargin{i*2});
            if(isnumeric(Arg))
                UClass=unique(Classes(A));
                if(Arg==Inf)
                    Arg=length(UClass);
                end
                CArg=UClass(Arg);
                B=~gIsEmpty(regexp(Classes,CArg));
            else
                GoodOnes=cellfun('isclass', Classes, 'char');
                B=zeros(size(Classes))>1;
                if(Arg(1)=='~')
                    C=gIsEmpty(regexp(Classes(GoodOnes),Arg(2:end),'ignorecase'));
                else
                    C=~gIsEmpty(regexp(Classes(GoodOnes),Arg,'ignorecase'));
                end
                B(GoodOnes)=C;
            end
        case LCats
            FieldName=Cats(strcmpi(CurCat,Cats));
            Classes={ShortInfosArr.(FieldName{1})};
            %             B=strcmp(Classes,varargin{i*2});
            if(isnumeric(Arg))
                if(all(cell2mat(cellfunc(@isnumeric,Classes))))
                    Classes=[ShortInfosArr.(FieldName{1})];
                    B=(Classes==Arg);
                else
                    UClass=unique(Classes(A));
                    if(isempty(UClass))
                        B=zeros(size(Classes))>1;
                    else
                        CArg=UClass(Arg);
                        B=~gIsEmpty(regexp(Classes,CArg));
                    end
                end
            else
                B=strcmpi(Classes,Arg);
            end
        otherwise
            disp('Error getSpecificInfos');
            error('Error getSpecificInfos');
    end
    A=A & B;
end
[Tmp Order]=sort({ShortInfosArr(A).SeriesDateTime});
F=find(A);
OutList=F(Order);
Out=ShortInfosArr(OutList);
% disp(['Found ' num2str(length(F)) ' scans']);