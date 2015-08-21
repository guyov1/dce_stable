function [AIF,std_AIF]=calc_AIF_from_curves(concentration_4D,inds4aif_mat,method)

    
    curves4aif_num=size(inds4aif_mat,1);
    curves4aif_mat=zeros(curves4aif_num,size(concentration_4D,4));
            
    for ii=1:curves4aif_num
        row=inds4aif_mat(ii,1);
        col=inds4aif_mat(ii,2);
        slice=inds4aif_mat(ii,3);
        curves4aif_mat(ii,:)=squeeze(concentration_4D(row,col,slice,:));
    end
    
    if strcmp(method,'average') && ~isempty(curves4aif_mat)
        AIF=mean(curves4aif_mat,1);
        std_AIF=std(curves4aif_mat,1,1);
    elseif isempty(curves4aif_mat)
        AIF=zeros(1,size(curves4aif_mat,2));
        std_AIF=zeros(1,size(curves4aif_mat,2));
    end
end