BaseP='C:\DATA\DSCTemporalRes_Phantom\';
D=dir(BaseP);
D=D([D.isdir]);
D=D(3:end);
DN={D.name}';
%%
for d=1:numel(DN)
    CurP=[BaseP DN{d} filesep];
%     gDicom2Nifti(CurP,[BaseP DN{d}(1:9) '.nii']);
    A=loadniidata([BaseP DN{d}(1:9) '_0001.nii']);
    nSlices(d)=size(A,3);
    D2=dir(CurP);
    CurFullInfo=dicominfo([CurP D2(3).name]);
    FullInfos(d)=CurFullInfo;
    TimeBetweenDCEVolsPerSlice=double(min(CurFullInfo.AcquisitionMatrix(CurFullInfo.AcquisitionMatrix>0))*CurFullInfo.RepetitionTime)/(1000*max(CurFullInfo.EchoTrainLength,1))*CurFullInfo.PercentPhaseFieldOfView/100.;
    TimeBetweenDCEVols(d)=TimeBetweenDCEVolsPerSlice*nSlices(d);
end
%%
for d=1:numel(DN)
    TrueTimes(d)=str2num(DN{d}(20:end-1));
end
OtherMethod=([FullInfos.Private_0019_105a]./[FullInfos.NumberOfTemporalPositions])/1e6;
%%
[TrueTimes; TimeBetweenDCEVols; OtherMethod]