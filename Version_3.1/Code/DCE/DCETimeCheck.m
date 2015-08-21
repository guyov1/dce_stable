% PT{1}='/data/Gilad/Brain-Tumor-Project/ACRIN_0625/1-Voll_Sara_Case38/1-Baseline_13_8_07/Study20070813_114944/Series20070813_124218 DCE main';
% PT{2}='/data/Gilad/Brain-Tumor-Project/ACRIN_0625/1-Voll_Sara_Case38/3-W8_28_10_07/Study20071028_100029/Series20071028_110337 DCE Main_FA23_partial-sequence';
% PT{3}='/data/Gilad/Brain-Tumor-Project/ACRIN_0625/5-Kalfus_Ronen_Case50/ROK_Baseline-11_9_07/study/Series20070911_131721 DCE main';
% PT{4}='/data/Gilad/Brain-Tumor-Project/ACRIN_0625/5-Kalfus_Ronen_Case50/ROK_W2-30_9_07/Study20070930_110748/Series20070930_115303 DCE main';
PT{1}='/data/Gilad/Brain-Tumor-Project/Health-Control/Full-Protocol/1SHAHAF_TOMRY/Study20081116_103251/SHTO_Se18_DCE main-after 24 sec';
%%
for i=1:numel(PT)
    D=dir([PT{i} '/*.dcm']);
    DCMFile{i}=[PT{i} filesep D(1).name];
    CurFullInfoC{i}=dicominfo(DCMFile{i});
    CurFullInfo=CurFullInfoC{i};
    CurFullInfos(i)=CurFullInfo;
    TimeBetweenDCEVolsPerSlice(i)=double(min(CurFullInfo.AcquisitionMatrix(CurFullInfo.AcquisitionMatrix>0))*CurFullInfo.RepetitionTime)/(1000*max(CurFullInfo.EchoTrainLength,1))*CurFullInfo.PercentPhaseFieldOfView/100.;
    
    OutP{i}=['/data/home/gilad/Desktop/DCE/gilad/Database/DCEOut/' ToShortName(CurFullInfo.PatientName.FamilyName) '_' CurFullInfo.SeriesDate filesep];
%     MeanVol=loadniidata([OutP{i} 'DCEMean.nii']);
%     SzC{i}=size(MeanVol);
%     TimeBetweenDCEVols(i)=TimeBetweenDCEVolsPerSlice(i)*SzC{i}(3);
        TimeBetweenDCEVols(i)=TimeBetweenDCEVolsPerSlice(i)*12;
end
