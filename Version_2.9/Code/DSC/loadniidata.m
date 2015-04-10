% function [Out Header A]=loadniidata(FN)
function [Out Header A]=loadniidata(FN, adsasd)
% if exist('load_untouch_nii','file')
%     A=load_untouch_nii(FN);
% else
%     A=load_nii(FN);
% end
try
    A=load_untouch_nii(FN);
catch
    A=load_nii(FN);
end
Out=double(A.img);
Header= rmfield(A, 'img');