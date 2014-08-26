function [Ct Keep_in_Mask] =Intens2concentration_voxel(intensity_curve,baseline_edges,TE)
% transform intensity curve to concentration curve ( C(t) ), up to a
% constant.
Keep_in_Mask=1;

s0=mean(intensity_curve(baseline_edges(1):baseline_edges(end)));
Ct=-log(intensity_curve/s0)/TE;
Ct=Ct(:);

% Correcting cases with C(t)= "Inf":

inf_inds_vec=find(Ct==Inf);

if ~isempty(inf_inds_vec)
    
    if length(inf_inds_vec)>1
        inf_inds_diff=inf_inds_vec(2:end)-inf_inds_vec(1:end-1);
    else
        inf_inds_diff=[];
    end
    
    if any(inf_inds_diff==1) %Throw C(t) with 2 or more consecutive "Inf":
        Keep_in_Mask=0;
    else % only singular Inf points - interpolate
        for ii=1:length(inf_inds_vec)
            inf_ind=inf_inds_vec(ii);
            if inf_ind==1 %first index
                Ct(inf_ind)=Ct(inf_ind+1);
            elseif inf_ind==length(Ct) % last index
                Ct(inf_ind)=Ct(inf_ind-1);
            else % all regular cases:
                Ct(inf_ind)=0.5* (Ct(inf_ind-1)+Ct(inf_ind+1));
            end
        end % for
    end
end


end % function


