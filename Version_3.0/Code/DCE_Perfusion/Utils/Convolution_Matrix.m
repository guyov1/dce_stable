function [ Conv_Matrix ] = Convolution_Matrix( min_interval, AIF )
%Convolution_Matrix Creates convolution matrix of AIF

% Make sure AIF is a vector
if ( size(AIF,2) ~= 1)
    AIF = AIF';
end

% ---------------------------------------------
% Piece-wise constant
% ---------------------------------------------

% num_time_stamps                 = length(AIF);
% [X, Y]                          = meshgrid(1:num_time_stamps);
% Out                             = ((Y-X)+1).*(X<=Y);
% Out(Out<=0)                     = num_time_stamps + 1;
% AIF_NaN_Pad                     = [AIF ; NaN];
% Conv_Matrix                     = min_interval*AIF_NaN_Pad(Out);
% Conv_Matrix(isnan(Conv_Matrix)) = 0;

% ---------------------------------------------
% Linear change - elegant implementation
% ---------------------------------------------
num_time_stamps                   = length(AIF);
[X, Y]                            = meshgrid(1:num_time_stamps);

Out                               = ((Y-X)+1).*(X<=Y);
% Out_plus_1                        = Out - 1;
% Out_minus_1                       = Out - 2;
Out_plus_1                        = Out + 1;
Out_minus_1                       = Out - 1;

Out(Out<=0)                       = num_time_stamps + 1;
Out_plus_1(Out_plus_1<=0)         = num_time_stamps + 1;
Out_minus_1(Out_minus_1<=0)       = num_time_stamps + 1;

% AIF_shifted = [AIF(2:end) ; 0];
% AIF_shifted = [0 ; AIF(1:end-1)];
% AIF_NaN_Pad                       = [AIF_shifted ; NaN];
AIF_NaN_Pad                       = [AIF ; NaN];

Term                              = AIF_NaN_Pad(Out);
NaN_Matrix                        = Term ./ Term ; % Used to mask unwanted indices
Term(isnan(Term))                 = 0;
Term_plus_1                       = AIF_NaN_Pad(Out_plus_1)  .* NaN_Matrix;
Term_plus_1(isnan(Term_plus_1))   = 0;
Term_minus_1                      = AIF_NaN_Pad(Out_minus_1) .* NaN_Matrix;
Term_minus_1(isnan(Term_minus_1)) = 0;
% Conv_Matrix                       = (1/6)*min_interval*( Term_minus_1 + 4*Term_plus_1 + Term);
Conv_Matrix                       = (1/6)*min_interval*( Term_minus_1 + 4*Term + Term_plus_1 );

% ---------------------------------------------
% Linear change - non-elegant implementation
% ---------------------------------------------
% num_time_stamps                 = length(AIF);
% Conv_Matrix = zeros(num_time_stamps);
% for i = 1 : num_time_stamps
%     for j = 1 : i
% for i = 0 : num_time_stamps-1
%     for j = 0 : i
%         First_valid = (i-j-1 >= 1) && (i-j-1 <= num_time_stamps);
%         Sec_valid   = (i-j   >= 1) && (i-j   <= num_time_stamps);
%         Third_valid = (i-j+1 >= 1) && (i-j+1 <= num_time_stamps);
%         
%         if (First_valid && Sec_valid && Third_valid)
%             Conv_Matrix(i,j) = min_interval* ( AIF(i-j-1) + 4*AIF(i-j) + AIF(i-j+1) ) / 6;
%         else if (~First_valid && Sec_valid && Third_valid)
%                 Conv_Matrix(i,j) = min_interval* ( 4*AIF(i-j) + AIF(i-j+1) ) / 6;
%             else if (~First_valid && ~Sec_valid && Third_valid)
%                     Conv_Matrix(i,j) = min_interval* ( AIF(i-j+1) ) / 6;
%                 else if (~First_valid && ~Sec_valid && ~Third_valid)
%                         Conv_Matrix(i,j) = 0;
%                     else if (First_valid && ~Sec_valid && ~Third_valid)
%                             Conv_Matrix(i,j) = min_interval* ( AIF(i-j-1) ) / 6;
%                         else if(First_valid && Sec_valid && ~Third_valid)
%                                 Conv_Matrix(i,j) = min_interval* ( AIF(i-j-1) + 4*AIF(i-j)) / 6;
%                             else if(First_valid && ~Sec_valid && Third_valid)
%                                     Conv_Matrix(i,j) = min_interval* (AIF(i-j-1) + AIF(i-j+1)) / 6;
%                                 else if(~First_valid && Sec_valid && ~Third_valid)
%                                         Conv_Matrix(i,j) = min_interval* ( 4*AIF(i-j) ) / 6;
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     
% end

end

