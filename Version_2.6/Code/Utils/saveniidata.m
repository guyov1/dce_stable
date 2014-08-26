function saveniidata(NII,FN)
% function saveniidata(NII,FN)
% if exist('save_untouch_nii','file')
%     save_untouch_nii(NII,FN);
% else
%     save_nii(NII,FN);
% end
try
    save_nii(NII,FN);
catch
    save_untouch_nii(NII,FN);
end