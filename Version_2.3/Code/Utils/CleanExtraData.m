% BaseP='C:\STRDCE\John\Database\DCEOut\';
BaseP='C:\DCE\John\Database\DCEOut\';
BaseP='C:\DCET3\gilad\Database\DCEOut\';
D=dir(BaseP);
D=D(3:end);
%%
for d=1:numel(D)
    disp(d);
    disp(D(d).name);
    CurP=[BaseP D(d).name filesep 'DCEMainCoreged\'];
    DD=dir([CurP '*.nii']);
    for dd=2:numel(DD)
        delete([CurP DD(dd).name]);
    end
    CurP=[BaseP D(d).name filesep 'DCEMainNii\'];
    DD=dir([CurP '*.nii']);
    for dd=2:numel(DD)
        delete([CurP DD(dd).name]);
    end
end