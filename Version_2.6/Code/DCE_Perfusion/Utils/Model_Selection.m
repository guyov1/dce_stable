function [ ChosenByAIC ] = Model_Selection( nDataPoints, RMSmaps, NParams, Data_Weight, AIC_Correction)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Input - 
%         nDataPoints - Number of volumes (originally nVols)
%         RMSmaps     - RMS maps - fitting to data which replaces
%                       a likelihood function
%         NParams     - Number of parameters in each method

% WorkingP='C:\STRDCE\John\Database\DCEOut\RuYe_20091224\';
% DataP=[WorkingP 'AutoArtBAT' filesep];
% a=load([WorkingP 'PKM3D.mat']);

RSSmaps = RMSmaps.^2 * nDataPoints; 
%RSSmaps = RMSmaps;

%%
replicated_NParams        = zeros(size(RSSmaps));
replicated_nDataPoints    = zeros(size(RSSmaps));
replicated_nDataPoints(:) = nDataPoints ;
for i = 1 : length(NParams)
    replicated_NParams(:,:,:,i) = NParams(i);
end

% MAD=median(abs(diff(CTC2DBigGood,[],2)),2);
% SSig=(MAD/0.67).^2;.*repmat(SSig,[1 4])
AIC          =  2 * replicated_NParams + Data_Weight * nDataPoints * log( RSSmaps );
%AIC          = ( 0 * log(nDataPoints) + 1 * 2 ) * repmat( NParams, [ size(RMSmaps, 1) 1] ) + nDataPoints * 0.1 * log( RMSmaps );
if AIC_Correction
    Correction   = 2 * replicated_NParams .* (replicated_NParams + 1) ./ ( replicated_nDataPoints - replicated_NParams - 1);
else
    Correction   = zeros(size(RSSmaps));
end

% Corrected Akaike
AICc         = AIC + Correction;

[~, ChosenByAIC] = min(AICc,[],4);

%ChosenParamMap = ParamMaps(:, :, :, ChosenByAIC);

% Tmp3D(MskX) = ChosenByAIC;
% 
% figure(1122);
% clf;
% montage(permute(mritransform(Tmp3D(GoodRows,GoodCols,GoodSlices)+1),[1 2 4 3]),[0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0]);
% colorbar('YTick',1:5,'YTickLabel',{'No computation','Empty model','Only plasma','Patlak','Tofts'});
% set(gcf,'Position',figposition([0 0 100 100]));
% title('AICc');
% saveas(1122,[WorkingP VPstr 'AICc.png']);
% saveas(1122,[WorkingP VPstr 'AICc.fig']);
% close(1122);
% 
% AddToLog(WorkingP,['ye' VPstr '_t21AICc'],[VPstr ' AICc.'],[VPstr 'AICc.png']);
% 
% %%
% Raw2Nii(Tmp3D,[PKOutP 'AICc' '.nii'],'float32', MeanFN);

end

