function  Rt_4D = SVD_solve_4D(AIF,Ct_4D,Mask,deltaT,Rt_type,deconv_methods,handles)

% This function solve (Ab=c) using SVD for b, (as shown in the paper of
% Ostergard 1996), for the whole 4D data (all voxels that passed the mask)
%   A is matrix built from AIF,
%   Ct=c=concentration curve at a voxel, method='constant' or 'linear' (assumption on the
%       difference between 2 time points),
%   b=Rt (the desired solution).
%
[N_rows,N_cols,N_slices,N_time_points]=size(Ct_4D);

% initialize:
if deconv_methods.sSVD.en
    display('sSVD');
    Rt_4D_sSVD=zeros(N_rows,N_cols,N_slices,N_time_points);
end
if deconv_methods.cSVD.en
    display('cSVD');
    Rt_4D_cSVD=zeros(N_rows,N_cols,N_slices,N_time_points);
end
if deconv_methods.oSVD.en
    display('oSVD');
    Rt_4D_oSVD=zeros(N_rows,N_cols,N_slices,N_time_points);
end

% Make the deconvolution for each voxel:
status_string=cell(N_slices,1);
data_voxels_row_col_slice=handles.data_voxels_row_col_slice;
slice_curr=0;
for ii=1:size(data_voxels_row_col_slice,1)
    row=data_voxels_row_col_slice(ii,1);
    col=data_voxels_row_col_slice(ii,2);
    slice=data_voxels_row_col_slice(ii,3);
    if slice~=slice_curr
        display(['Deconvolving Slice #',num2str(slice)]);
    end
    Ct=Ct_4D(row,col,slice,:);
    [Rt]=SVD_solve_voxel(AIF,Ct(:),deltaT,Rt_type,deconv_methods);
    if deconv_methods.sSVD.en
        Rt_4D_sSVD(row,col,slice,:)=Rt.sSVD;
    end
    if deconv_methods.cSVD.en
        Rt_4D_cSVD(row,col,slice,:)=Rt.cSVD(1:N_time_points);
    end
    if deconv_methods.oSVD.en
        Rt_4D_oSVD(row,col,slice,:)=Rt.oSVD;
    end
    if (ii<size(data_voxels_row_col_slice,1) && data_voxels_row_col_slice(ii+1,3)~=slice) || ii==size(data_voxels_row_col_slice,1)
        status_string{slice}=['Slice ',num2str(slice),' ......... Done'];
        set(handles.deconvolution_status_textbox,'String',status_string);
    end
    slice_curr=slice;
end

% for slice=1:N_slices
%     display(['Deconvolving Slice #',num2str(slice)]);
%     for row=1:N_rows
%         parfor col=1:N_cols
%             if Mask(row,col,slice)==1
%                 Ct=Ct_4D(row,col,slice,:);
%                 [Rt]=SVD_solve_voxel(AIF,Ct(:),deltaT,Rt_type,deconv_methods);
%                 if deconv_methods.sSVD.en
%                     Rt_4D_sSVD(row,col,slice,:)=Rt.sSVD;
%                 end
%                 if deconv_methods.cSVD.en
%                     Rt_4D_cSVD(row,col,slice,:)=Rt.cSVD(1:N_time_points);
%                 end
%                 if deconv_methods.oSVD.en
%                     Rt_4D_oSVD(row,col,slice,:)=Rt.oSVD;
%                 end
%             end
%         end
%     end
%     status_string{slice}=['Slice ',num2str(slice),' ......... Done'];
%     set(handles.deconvolution_status_textbox,'String',status_string);
% end

% arrange the output:
if deconv_methods.sSVD.en
    Rt_4D.sSVD=Rt_4D_sSVD;
end
if deconv_methods.cSVD.en
    Rt_4D.cSVD=Rt_4D_cSVD;
end
if deconv_methods.oSVD.en
    Rt_4D.oSVD=Rt_4D_oSVD;
end
display('Done!')
