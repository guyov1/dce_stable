function [Flow,CBV,MTT,K1,K2] = calc_maps(concentration4D,AIF,Mask,deconv_methods,handles)

% Input:  concentration4D: c(t) for each voxel (that passed the mask)
%         AIF: was chosen manually or automotically before
%         Mask: equals '1' if voxel is relevant, '0' if doesn't
%
% Output: Flow: Ft for each voxel that is relevant
%
%   more parameters to come...

Rt_type='linear'; %(assumption that R(t) between 2 time points vary linearly)
deltaT=2; % [sec]. Should be taken from the data.

% calc R(t) (for all voxels passed filters):
% regularize (if user wanted)
% if deconv_methods.reg || deconv_methods.reg_SVD
%    % regularize the matrix equation... 
% end

data_voxels_indices=handles.data_voxels_indices;

%% calculating CBV
% Use eq.3 from the paper of Ostergaard-1996:
% CBV=integral(C_VOI(t)) / integral(AIF(t))
% CBV=zeros(size(concentration4D));
% CBV=sum(concentration4D,4)/sum(AIF);
CBV_init=trapz(concentration4D,4)/trapz(AIF);
CBV.corr=CBV_init;
CBV.no_corr=CBV_init;
display('Calculated CBV map');
% Inserting correction for contrast agent extravasation, according to the
% paper of Weisskoff 2006
% The defualt CBV will be the one WITH the correction.
% CBV_no_corr is the CBV WITHOUT the correction.
Enhancment_ratio_3D= (handles.Mask) .* (mean(handles.data4D_init(:,:,:,end-9:end),4)) ./ (handles.data_average_bl);
if handles.permeability_correction
    K1=zeros(size(CBV_init));
    K2=zeros(size(CBV_init));
    for slice=1:size(concentration4D,3)
        % Estimate noise:
        data3D_slice=handles.data4D_init(:,:,slice,:);
        conc3D_slice=concentration4D(:,:,slice,:);
%         data_var_around_zero=mean(data3D_slice.^2,4);
%         NonBrainVoxels=find(handles.BrainMask(:,:,slice)==0);
%         ZeroVoxels=find(data_var_around_zero==0);
%         NonZeroNoiseInds=setdiff(NonBrainVoxels,ZeroVoxels);
%         noise_std_slice=sqrt(mean(data_var_around_zero(NonZeroNoiseInds)));
%         [noise_voxels_rows noise_voxels_cols noise_voxels_slices]=ind2sub(size(CBV),);
        
        % We use our Mask (from pre-processing stage) as mask #1.
        Mask1=handles.Mask(:,:,slice);
        
        %calc mask #2 - of non-enhancing voxels
        average_bl_slice=handles.data_average_bl(:,:,slice);
%         std_over_average_bl_slice=average_bl_slice+std(data3D_slice(:,:,1,handles.first_bl_sample:handles.last_bl_sample),1,4);
%         std_below_average_bl_slice=average_bl_slice-5*std(data3D_slice(:,:,1,handles.first_bl_sample:handles.last_bl_sample),1,4);
        data_steady_state_avg=mean(data3D_slice(:,:,1,end-9:end),4);
%         NonEnhancingMask=Mask1.*(data_steady_state_avg>std_below_average_bl_slice);
        NonEnhancingMask=(data_steady_state_avg>0.9*average_bl_slice);
%         NonEnhancingMask=Mask1;
        PermCorrMask=NonEnhancingMask.*Mask1;
        [Mask_rows Mask_cols]=ind2sub(size(Mask1),find(Mask1>0));
        
        %calc 2 function of time - the average c(t), and its cumulative trapezoidal
        %sum
        conc3D_slice_masked=conc3D_slice .* repmat(PermCorrMask,[1 1 1 size(conc3D_slice,4)]);
        conc_mean_t=squeeze(sum(sum(conc3D_slice_masked))) ./ sum(sum(PermCorrMask));
        conc_mean_integral_t=2*cumtrapz(conc_mean_t); %calculate intergral using trapezoidal method. deltaT=2[sec]
        
        %find K1,K2 paramters by linear LS fitting:
        K1_slice=zeros(size(Mask1));
        K2_slice=zeros(size(Mask1));
        for ind=1:length(Mask_rows)
           row=Mask_rows(ind);
           col=Mask_cols(ind);
           Ct=squeeze(conc3D_slice(row,col,1,:));
           Kvec = inv( [conc_mean_t';conc_mean_integral_t']*[conc_mean_t conc_mean_integral_t])*[conc_mean_t';conc_mean_integral_t']*Ct;
           K1_slice(row,col)=Kvec(1); 
           K2_slice(row,col)=-Kvec(2);
        end
        K1(:,:,slice)=K1_slice;
        K2(:,:,slice)=K2_slice;
        %add the correction to the CBV (for each voxel):
        CBV.corr(:,:,slice)=CBV_init(:,:,slice)+K2_slice*trapz(conc_mean_integral_t);
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculating Flow and MTT- using SVD methods:
% MTT = CBV/Flow
if any([deconv_methods.sSVD.en deconv_methods.cSVD.en deconv_methods.oSVD.en ])
    display('Start calculating CBF using methods:');
    Rt_4D = SVD_solve_4D(AIF,concentration4D,Mask,deltaT,Rt_type,deconv_methods,handles);
    if deconv_methods.sSVD.en
        Flow.sSVD=max(Rt_4D.sSVD,[],4);
        MTT.sSVD=zeros(size(CBV_init));
        MTT.sSVD(data_voxels_indices)=CBV.corr(data_voxels_indices) ./ Flow.sSVD(data_voxels_indices);
    end
    if deconv_methods.cSVD.en    
        Flow.cSVD=max(Rt_4D.cSVD,[],4);
        MTT.cSVD=zeros(size(CBV_init));
        MTT.cSVD(data_voxels_indices)=CBV.corr(data_voxels_indices) ./ Flow.cSVD(data_voxels_indices);
    end
    if deconv_methods.oSVD.en
        Flow.oSVD=max(Rt_4D.oSVD,[],4);
        MTT.oSVD=zeros(size(CBV_init));
        MTT.oSVD(data_voxels_indices)=CBV.corr(data_voxels_indices) ./ Flow.oSVD(data_voxels_indices);
    end
end

