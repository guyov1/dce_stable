function [ OutName ] = NormalizeNii( Path, FileName, RefNii, Threshold, Threshold_val )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%RefNii        = DCEFNs{1};
% Threshold     = true;
% Threshold_val = 150;

NiiMap = loadniidata([Path filesep FileName '.nii']);

if Threshold
    max_val                                          = Threshold_val;
    NiiMap_Thresholded                               = NiiMap;
    NiiMap_Thresholded(NiiMap_Thresholded > max_val) = max_val;
    NiiMap_Norm_0_1                                  = NiiMap_Thresholded ./ max(max(max(NiiMap_Thresholded)));
    suffix                                           = ['_Thresholded_' num2str(Threshold_val) '_Normalized_0_1.nii'];
else
    NiiMap_Norm_0_1                                  = NiiMap ./ max(max(max(NiiMap)));
    suffix = '_Normalized_0_1.nii';
end


OutName = [Path filesep FileName suffix];

Raw2Nii(NiiMap_Norm_0_1, OutName, 'float32', RefNii);

end


