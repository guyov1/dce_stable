function Out=bwfillHoles3Dby2D(Out)
for i=1:size(Out,3)
    Out(:,:,i)=bwfill(squeeze(Out(:,:,i)),'holes');
end