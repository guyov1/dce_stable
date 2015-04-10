function [Ct] =voxel_intens2concentration(intensity_curve,baseline_edges,TE)
% transform intensity curve to concentration curve ( C(t) ), up to a
% constant.
s0=mean(intensity_curve(baseline_edges(1):baseline_edges(end)));
Ct=-log(intensity_curve/s0)/TE;

end


