function [Rval,  Pval] = CorrNiiMaps( Path_Map1, Path_Map2, SliceNum, MinRange, MaxRange, Mask1, Mask2, Mask3, Mask4, Use_Mask)


NiiMap1   = loadniidata(Path_Map1);
NiiMap2   = loadniidata(Path_Map2);

NiiMap1   = NiiMap1(:, :, SliceNum);
NiiMap2   = NiiMap2(:, :, SliceNum);

Mask = uint8(Mask1) | uint8(Mask2) | uint8(Mask3) | uint8(Mask4);

NiiMap1_Masked = zeros(size(NiiMap1));
NiiMap2_Masked = zeros(size(NiiMap2));

if Use_Mask
    NiiMap1_Masked(Mask(:, :, SliceNum)==1)   = NiiMap1(Mask(:, :, SliceNum)==1);
    NiiMap2_Masked(Mask(:, :, SliceNum)==1)   = NiiMap2(Mask(:, :, SliceNum)==1);
    
    figure;
    subplot(2,2,1);
    imshow(NiiMap1);
    title(['DCE Map. Slice: ' num2str(SliceNum)]);
    subplot(2,2,2);
    imshow(NiiMap1_Masked,'Colormap',jet(255));
    subplot(2,2,3);
    imshow(NiiMap2);
    title(['DSC Map. Slice: ' num2str(SliceNum)]);
    subplot(2,2,4);
    imshow(NiiMap2_Masked,'Colormap',jet(255));
else
    NiiMap1_Masked = NiiMap1;
    NiiMap2_Masked = NiiMap2;
end

maxVal1   = max(max(max(NiiMap1_Masked)));
maxVal2   = max(max(max(NiiMap2_Masked)));
minVal1   = min(min(min(NiiMap1_Masked)));
minVal2   = min(min(min(NiiMap2_Masked)));

idx1      = find( NiiMap1_Masked > (MinRange*minVal1) & NiiMap1_Masked < (MaxRange*maxVal1) );
idx2      = find( NiiMap2_Masked > (MinRange*minVal2) & NiiMap2_Masked < (MaxRange*maxVal2) );
final_idx = intersect(idx1, idx2);


vector1   = NiiMap1_Masked( final_idx);
vector2   = NiiMap2_Masked( final_idx);

%[R, PValue] = corrplot([vector1  vector2], 'testR','on');
[Rval, Pval] = corrplot([vector1  vector2], 'testR','on');

%SliceNum = 4 ; Path1 ='\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\Perfusion_DCE\coregFlow\DSC2DCE\Flow_Larsson_Relative_WM_30_6_brain_Thresholded_200_Normalized_0_1.nii' ; Path2 ='\\fmri-t9\users\Moran\DCE\HTR_STROKE\01_REMEZ_YECHEZKEL\Study20140615_114415\Perfusion_DCE\coregFlow\DSC2DCE\rdsc_oCBFlr_Thresholded_200_Normalized_0_1.nii' ;CorrNiiMaps(Path1, Path2,4);

end



