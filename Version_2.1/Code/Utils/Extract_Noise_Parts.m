function [ noise_matrix ] = Extract_Noise_Parts( in_image, brain_mask )
%Extract_Noise_Parts  - Extracts noise parts from brain image
%   The function tries to extract all the pixels that do not contain real
%   signal.
%   Input  - Image
%   Output - matrix of noise part

close all;

%% Vars initiatios

% Initiate general variables
rows           = size(in_image,1);
columns        = size(in_image,2);
middle_row_idx = round(rows/2);
middle_col_idx = round(columns/2);

% Num of bins for histogram
num_hist_bins = 55;

% Mask of all 0 value (will be ignored later)
Zero_Mask = find(in_image==0);

% Get the min and max values of image
Min_Val = min(min(in_image));
Max_Val = max(max(in_image));

% Zero extreme image values
in_image_zeroed = in_image;
in_image_zeroed(in_image_zeroed>0.95*(Max_Val - Min_Val)) = 0;
in_image_zeroed(in_image_zeroed>0.03*(Max_Val - Min_Val)) = 0;

% Create edge image
in_image_zeroed_edge = edge(in_image_zeroed);

% Display before and after zeroing effect
figure;
subplot(2,2,1);
imshow(mat2gray(in_image));
title('Image before zeroing extreme values');
subplot(2,2,2);
imshow(mat2gray(in_image_zeroed));
title('Image after zeroing extreme values');
subplot(2,2,3);
imshow(mat2gray(in_image_zeroed_edge));
title('Edge image after zeroing extreme values');

%% Get image boundaries (which contains data ~= 0)

% Get the x boundaries of the image
num_rows_to_check = 5;
all_min           = zeros(1,num_rows_to_check);
all_max           = zeros(1,num_rows_to_check);
iter_idx          = 1;

% Check 5 rows in the middle of the image for the average min and max
% indices from which there is a real signal value
for j = middle_row_idx-2:middle_row_idx+2
    % Go to the middle of the image and get all column indices with value>Min_Val
    % The assumption is that 0 value is not a real signal
    good_idx = find(image_to_work(j,:)>Min_Val);
    % Get the min and max indices of the image pixels
    min_idx = min(good_idx);
    max_idx = max(good_idx);
    
    all_min(iter_idx) = min_idx;
    all_max(iter_idx) = max_idx;
    
    iter_idx = iter_idx + 1;
end

% Average the min and the max columns
avg_min_col = round(mean(all_min));
avg_max_col = round(mean(all_max));

% Get the y boundaries of the image
num_cols_to_check = 5;
all_min = zeros(1,num_cols_to_check);
all_max = zeros(1,num_cols_to_check);
iter_idx = 1;

% Check 5 columns in the middle of the image for the average min and max
% indices from which there is a real signal value
for j = middle_col_idx-2:middle_col_idx+2
    % Go to the middle of the image and get all rows indices with value>Min_Val
    % The assumption is that 0 value is not a real signal
    good_idx = find(image_to_work(:,j)>Min_Val);
    min_idx = min(good_idx);
    max_idx = max(good_idx);
    % Get the min and max indices of the image pixels
    all_min(iter_idx) = min_idx;
    all_max(iter_idx) = max_idx;
    
    iter_idx = iter_idx + 1;
end

% Average the min and the max columns
avg_min_row = round(mean(all_min));
avg_max_row = round(mean(all_max));

%% Get the parts of the noise of the truncated image

% Get the edge map which will help us decide where is the noise boundaries
edge_map      = edge(image_to_work);
cols_to_check = round( (avg_min_col + avg_max_col) / 2);
total_cols    = 5;

possible_min_rows = zeros(1,total_cols);
possible_max_rows = zeros(1,total_cols);

iter_idx = 1;

for j = cols_to_check-2:cols_to_check+2
    possible_min_idx = find(edge_map(:,j)>0, 1 );
    possible_max_idx = find(edge_map(:,j)>0, 1, 'last' );
    
    possible_min_rows(iter_idx) = possible_min_idx;
    possible_max_rows(iter_idx) = possible_max_idx;
    
    iter_idx = iter_idx + 1;
end

% # of pixels to remove to be on the safe side
safe_side_margin = 6;

min_row_for_noise = round(mean(possible_min_rows)) - safe_side_margin;
max_row_for_noise = round(mean(possible_max_rows)) + safe_side_margin;


% Get truncated noise parts
trunc_noise_1 = image_to_work(avg_min_row:min_row_for_noise,avg_min_col:avg_max_col);
trunc_noise_2 = image_to_work(max_row_for_noise:avg_max_row,avg_min_col:avg_max_col);

total_noise = [ trunc_noise_1(:) ; trunc_noise_2(:) ];
total_noise = double(total_noise);


%% Display truncated image and noise parts

% Display truncated image
truncated_image = image_to_work(avg_min_row:avg_max_row,avg_min_col:avg_max_col);
figure;
subplot(3,2,1);
imshow(mat2gray(truncated_image));
title('Truncated Image');
subplot(3,2,2);
imshow(mat2gray(image_to_work));
title('Original Image');

% Display Trunacted noise
subplot(3,2,3);
imshow(mat2gray(trunc_noise_1));
title('Upper noise part');
subplot(3,2,4);
imshow(mat2gray(trunc_noise_2));
title('Lower noise part');

% Display histogram
subplot(3,2,5);
hist(total_noise,num_hist_bins);
title('Image Histogram');

check=double(image_to_work(230:250,40:218));check=check(:);

all_histogram = hist(total_noise,num_hist_bins);


%% Display all histogram next to each other

figure;
subplot(1,3,1);
bar(all_histogram);
title(' Noise Histogram');

noise_matrix = total_noise;

end




