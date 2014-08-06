function TMat=getTMatFromNii(FN)
[A B]=loadniidata(FN);
TMat=eye(4);
TMat(1,:)=B.hdr.hist.srow_x;
TMat(2,:)=B.hdr.hist.srow_y;
TMat(3,:)=B.hdr.hist.srow_z;