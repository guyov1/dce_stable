function Out=getComputerParams(ParamName, varargin)
Out=getComputerParamM(ParamName, varargin);
% % BaseP='/data/Niftis2/';
% BaseP='D:\users\gilad\Database\';
% switch lower(ParamName)
%     case 'projectnamepathdepth'
%         Out=4;
%     case 'basepath'
%         Out=BaseP;
%     case 'temppath'
%         Out=[BaseP 'Temp'];
%     case 'infosfn'
%         Out=[BaseP 'Infos' filesep 'Infos.mat'];
%     case 'remarkinfosfn'
%         Out=[BaseP 'Infos' filesep 'RemarkInfos.mat'];
%     case 'shortinfosfn'
%         Out=[BaseP 'Infos' filesep 'ShortInfos.mat'];
%     case 'veryshortinfosfn'
%         Out=[BaseP 'Infos' filesep 'T9_VeryShortInfos.mat'];
%     case 'name'
%         Out='T9';
%     case 'tpm'
%         Out={'/data/spm8/tpm/grey.nii';'/data/spm8/tpm/white.nii';'/data/spm8/tpm/csf.nii'};
%     case 'nniftispaths'
%         Out=1;
%     case 'niftispath'
%         NiftisPath{1}='/data/Niftis/';
% %         NiftisPath{2}='H:\Niftis\';
%         Out=NiftisPath{varargin{1}};
%     case 'basedir'
%         Out='/data/Gilad/Base';
%     otherwise
%         disp('Error getComputerParams');
%         error('Error getComputerParams');
% end