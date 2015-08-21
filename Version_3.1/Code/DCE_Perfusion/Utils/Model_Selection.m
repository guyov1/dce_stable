function [ ChosenByAIC ] = Model_Selection( nDataPoints, RMSmaps, NParams, Data_Weight, AIC_Correction)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Input - 
%         nDataPoints - Number of volumes (originally nVols)
%         RMSmaps     - RMS maps - fitting to data which replaces
%                       a likelihood function
%         NParams     - Number of parameters in each method


% Convert RMS to RSS
RSSmaps                   = RMSmaps.^2 * nDataPoints; 

replicated_NParams        = zeros(size(RSSmaps));
replicated_nDataPoints    = zeros(size(RSSmaps));
replicated_nDataPoints(:) = nDataPoints ;
for i = 1 : length(NParams)
    replicated_NParams(:,:,:,i) = NParams(i);
end

% MAD  = median(abs(diff(CTC2DBigGood,[],2)),2);
% SSig = (MAD/0.67).^2;.*repmat(SSig,[1 4])
% AIC         = ( 0 * log(nDataPoints) + 1 * 2 ) * repmat( NParams, [ size(RMSmaps, 1) 1] ) + nDataPoints * 0.1 * log( RMSmaps );
AIC          =  2 * replicated_NParams + Data_Weight * nDataPoints * log( RSSmaps );

if AIC_Correction
    Correction   = 2 * replicated_NParams .* (replicated_NParams + 1) ./ ( replicated_nDataPoints - replicated_NParams - 1);
else
    Correction   = zeros(size(RSSmaps));
end

% Corrected Akaike
AICc         = AIC + Correction;

[~, ChosenByAIC] = min(AICc,[],4);

%ChosenParamMap = ParamMaps(:, :, :, ChosenByAIC);

end

