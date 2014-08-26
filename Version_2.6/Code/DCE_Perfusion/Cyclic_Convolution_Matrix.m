function [ Conv_Matrix ] = Cyclic_Convolution_Matrix( min_interval, AIF)
%Convolution_Matrix Creates cyclic convolution matrix of AIF
% The padding is to the size of M+N+1 (M - length of AIF, filter - length of filter kernel)


% ---------------------------------------------
% Piece-wise constant
% ---------------------------------------------
% nTimePoints_padded              = 2 * length(AIF) + 1;
% AIF_Padded                      = [AIF ;zeros(length(AIF) + 1,1)];
% 
% [X, Y]                          = meshgrid(1:nTimePoints_padded);
% Out                             = ((Y-X)+1).*(X<=Y);
% Out(Out<=0)                     = nTimePoints_padded + 1;
% 
% AIF_NaN_Pad                     = [AIF_Padded ; NaN];
% Conv_Matrix                     = min_interval*AIF_NaN_Pad(Out);
% Conv_Matrix(isnan(Conv_Matrix)) = 0;
% 
% % Make the matrix circular
% L = nTimePoints_padded;
% for i = 1:L
%     for j = 1:L
%         if j>i
%             Conv_Matrix(i,j) = Conv_Matrix(L+i-j+1,1);
%         end
%     end
% end



% ---------------------------------------------
% Linear change - elegant implementation
% ---------------------------------------------
nTimePoints_padded               = 2 * length(AIF) + 1;
AIF_Padded                       = [AIF ;zeros(length(AIF) + 1,1)];

[X, Y]                            = meshgrid(1:nTimePoints_padded);

Out                               = ((Y-X)+1).*(X<=Y);
% Out_plus_1                        = Out - 1;
% Out_minus_1                       = Out - 2;
Out_plus_1                        = Out + 1;
Out_minus_1                       = Out - 1;

Out(Out<=0)                       = nTimePoints_padded + 1;
Out_plus_1(Out_plus_1<=0)         = nTimePoints_padded + 1;
Out_minus_1(Out_minus_1<=0)       = nTimePoints_padded + 1;

% AIF_shifted = [AIF(2:end) ; 0];
% AIF_shifted = [0 ; AIF(1:end-1)];
% AIF_NaN_Pad                       = [AIF_shifted ; NaN];
AIF_NaN_Pad                       = [AIF_Padded ; NaN];

Term                              = AIF_NaN_Pad(Out);
NaN_Matrix                        = Term ./ Term ; % Used to mask unwanted indices
Term(isnan(Term))                 = 0;
Term_plus_1                       = AIF_NaN_Pad(Out_plus_1)  .* NaN_Matrix;
Term_plus_1(isnan(Term_plus_1))   = 0;
Term_minus_1                      = AIF_NaN_Pad(Out_minus_1) .* NaN_Matrix;
Term_minus_1(isnan(Term_minus_1)) = 0;
% Conv_Matrix                       = (1/6)*min_interval*( Term_minus_1 + 4*Term_plus_1 + Term);
Conv_Matrix                       = (1/6)*min_interval*( Term_minus_1 + 4*Term + Term_plus_1 );

% Make the matrix circular
L = nTimePoints_padded;
for i = 1:L
    for j = 1:L
        if j>i
            Conv_Matrix(i,j) = Conv_Matrix(L+i-j+1,1);
        end
    end
end

end

