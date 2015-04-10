function Out=GetParameterFromSiemensMrProt(Prot,str)
% Prot=char(Infos.A1_3_12_2_1107_5_2_19_45773_2014112512282752183025365_0_0_0.Private_0029_1020');
Out=-778;
I=strfind(Prot,str);
if(isempty(I)) return; end
tmp=regexp(Prot(I(1):I(1)+200),'=\W*(\d*)\W*','tokens');
Out=str2num(tmp{1}{1});