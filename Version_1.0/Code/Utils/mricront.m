function mricront(M,M2)
Raw2Nii(M,'/tmp/XX.nii','float32');
if(nargin==1)
    mricronx('/tmp/XX.nii');
    return;
end
Raw2Nii(M2,'/tmp/XX2.nii','float32');
mricronx({'/tmp/XX.nii','/tmp/XX2.nii'},[0 1]);