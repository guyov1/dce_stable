function [full_maps_dir_string]  = save_maps(CBF,CBV,MTT,K1,K2,Mask,deconv_methods,dest_dir)

CBF_in=CBF;
CBV_in=CBV;
MTT_in=MTT;
K1_in=K1;
K2_in=K2;

maps_dir_string=['maps_',datestr(now,'dd_mm_yyyy_HHMMSS')];
% curr_dir=handles.dest_dir;
res_dir=[dest_dir,filesep,'results'];
% cd(dest_dir);
if ~isdir(res_dir)
    mkdir(res_dir);
end
% cd 'results';
% results_dir=pwd;
full_maps_dir_string=[res_dir,filesep,maps_dir_string];
mkdir(full_maps_dir_string);
% cd(curr_dir);

params_vec={'CBF','CBV','MTT','K1','K2'};

% save CBV, K1, K2, Mask - which are not dependent on the SVD method:
% Also rotate back to a "nii" position
CBV_in_corr=CBV_in.corr;
CBV_in_no_corr=CBV_in.no_corr;
CBV_nii_corr=make_nii(permute(flipdim(CBV_in_corr,1),[2 1 3 4]));
save_nii(CBV_nii_corr,[full_maps_dir_string,filesep,'CBV']);
save([full_maps_dir_string,filesep,'CBV.mat'],'CBV_in_corr');
CBV_nii_no_corr=make_nii(permute(flipdim(CBV_in_no_corr,1),[2 1 3 4]));
save_nii(CBV_nii_no_corr,[full_maps_dir_string,filesep,'CBV_no_corr']);
save([full_maps_dir_string,filesep,'CBV_no_corr.mat'],'CBV_in_no_corr');

K1_nii=make_nii(permute(flipdim(K1_in,1),[2 1 3 4]));
K2_nii=make_nii(permute(flipdim(K2_in,1),[2 1 3 4]));
save_nii(K1_nii,[full_maps_dir_string,filesep,'K1']);
save_nii(K2_nii,[full_maps_dir_string,filesep,'K2']);
save([full_maps_dir_string,filesep,'K1.mat'],'K1_in');
save([full_maps_dir_string,filesep,'K2.mat'],'K2_in');

Mask_nii=make_nii(permute(flipdim(double(Mask),1),[2 1 3]));
save_nii(Mask_nii,[full_maps_dir_string,filesep,'Mask']);
save([full_maps_dir_string,filesep,'Mask.mat'],'Mask');

% save CBF,MTT - which depend on the SVD method:
if deconv_methods.sSVD.en
    CBF=CBF_in.sSVD;
    MTT=MTT_in.sSVD;
    %rotation is done in the "make_nii" line
    
    %saving:
    for ii=1:length(params_vec) 
        if ~( strcmp(params_vec{ii},'CBV') || strcmp(params_vec{ii},'K1') || strcmp(params_vec{ii},'K2') )
            param=eval(params_vec{ii});
            param_nii=make_nii(permute(flipdim(param,1),[2 1 3 4]));
            save_nii(param_nii,[full_maps_dir_string,filesep,params_vec{ii},'_sSVD']);
            save([full_maps_dir_string,filesep,params_vec{ii},'_sSVD.mat'],params_vec{ii});
        end
    end
end

if deconv_methods.cSVD.en
    CBF=CBF_in.cSVD;
    MTT=MTT_in.cSVD;
    %rotation is done in the "make_nii" line
    
    %saving:
    for ii=1:length(params_vec) 
        if ~( strcmp(params_vec{ii},'CBV') || strcmp(params_vec{ii},'K1') || strcmp(params_vec{ii},'K2') )
            param=eval(params_vec{ii});
            param_nii=make_nii(permute(flipdim(param,1),[2 1 3 4]));
            save_nii(param_nii,[full_maps_dir_string,filesep,params_vec{ii},'_cSVD']);
            save([full_maps_dir_string,filesep,params_vec{ii},'_cSVD.mat'],params_vec{ii});
        end
    end
end

if deconv_methods.oSVD.en
    CBF=CBF_in.oSVD;
    MTT=MTT_in.oSVD;
    %rotation is done in the "make_nii" line
    
    %saving:
    for ii=1:length(params_vec) 
        if ~( strcmp(params_vec{ii},'CBV') || strcmp(params_vec{ii},'K1') || strcmp(params_vec{ii},'K2') )
            param=eval(params_vec{ii});
            param_nii=make_nii(permute(flipdim(param,1),[2 1 3 4]));
            save_nii(param_nii,[full_maps_dir_string,filesep,params_vec{ii},'_oSVD']);
            save([full_maps_dir_string,filesep,params_vec{ii},'_oSVD.mat'],params_vec{ii});            
        end
    end
end
