function insertTMatToNii(FN,TMat)
[Q QQ Nii]=loadniidata(FN);
Nii.hdr.hist.sform_code=2;
Nii.hdr.hist.srow_x=TMat(1,:);
Nii.hdr.hist.srow_y=TMat(2,:);
Nii.hdr.hist.srow_z=TMat(3,:);
saveniidata(Nii,FN);