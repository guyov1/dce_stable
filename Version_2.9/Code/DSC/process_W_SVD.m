function [W_inv_out] = process_W_SVD(W_inv_in,th)

% W_in and W_out are diagonal matrices, made by SVD
% in the processing we usually eliminate small values in W, to "regularize"
% the solution of Ab=c.
W_inv_out=W_inv_in;

th_perc=th; % percentage of max element in W. Should be taken from GUI
max_element_in_W=1/min(min(diag(W_inv_in))); %max element in W is the min element in inv(W) (they are diagonal)
th_value=max_element_in_W*th_perc/100;
W_inv_out(find(W_inv_out>1/th_value))=0;


