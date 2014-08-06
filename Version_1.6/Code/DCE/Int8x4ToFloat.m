TmpFN=[WorkingP 'a.a'];
fid=fopen(TmpFN,'w')
fwrite(fid,CurFullInfo.Private_0019_105a','uint8')
fclose(fid);
fid=fopen(TmpFN);
A=fread(fid,1,'float32')
fclose(fid);