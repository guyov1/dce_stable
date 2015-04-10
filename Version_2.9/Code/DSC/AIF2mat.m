function [A] = AIF2mat(AIF,deltaT,method,params)
% This function build the AIF time-curve as a matrix.
% Input:    AIF
%           params - such as delta t, start point, end point...
%           method - 'constant' (assume concentraions are constant between
%                               time points)
%                     or:
%                    'linear' (assume concentraions vary linearly between
%                               time points)

N=length(AIF);
A=zeros(N,N);
if strcmp(method,'constant')
    for ii=1:N
        for jj=1:ii
            A(ii,jj)=AIF(ii-jj+1);
        end
    end
    
elseif strcmp(method,'linear')
    for ii=1:N
        for jj=1:ii
%             first_ind=(ii-jj)+1;  %in the paper it's ii-jj-1, but then the diagonal need AIF(-1), so I take AIF(0) - in matlab it's AIF(1)
%             % regular case of A(i,j) - not in boundaries:
%             if first_ind+2<=length(AIF) 
%                 A(ii,jj)=1/6 * (AIF(first_ind)+4*AIF(first_ind+1)+AIF(first_ind+2));
% %             boundary case #1 - ind2 out of bounds, ind1 in bounds 
%             elseif first_ind+1<=length(AIF)
%                 A(ii,jj)=1/5 * (AIF(first_ind)+4*AIF(first_ind+1));
% %             boundary case #2 - ind2 and ind1 out of bounds:                 
%             else
%                 A(ii,jj)=AIF(first_ind);
%             end
            main_ind=ii-jj+1;
            % regular case - main index is not in boundries, it has 2 non-zero neighbors
            if 1<main_ind && main_ind<length(AIF)
               A(ii,jj)= 1/6 * (AIF(main_ind-1)+4*AIF(main_ind)+AIF(main_ind+1));
            % main index is the first one. (it has only one neighbor)
            elseif main_ind==1
                A(ii,jj)= 1/5 * (4*AIF(main_ind)+AIF(main_ind+1));
            % main index is the last one. (it has only one neighbor)
            else
               A(ii,jj)= 1/5 * (AIF(main_ind-1)+4*AIF(main_ind));
            end
        end
    end
end

A=A*deltaT;