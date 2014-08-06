function [ Conv_Matrix ] = Cyclic_Convolution_Matrix( min_interval, AIF)
%Convolution_Matrix Creates cyclic convolution matrix of AIF
% The padding is to the size of M+N+1 (M - length of AIF, filter - length of filter kernel)


nTimePoints_padded              = 2 * length(AIF) + 1;
AIF_Padded                      = [AIF ;zeros(length(AIF) + 1,1)];

[X, Y]                          = meshgrid(1:nTimePoints_padded);
Out                             = ((Y-X)+1).*(X<=Y);
Out(Out<=0)                     = nTimePoints_padded + 1;

AIF_NaN_Pad                     = [AIF_Padded ; NaN];
Conv_Matrix                     = min_interval*AIF_NaN_Pad(Out);
Conv_Matrix(isnan(Conv_Matrix)) = 0;

end

