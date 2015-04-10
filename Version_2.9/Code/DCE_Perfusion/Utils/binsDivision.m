function [real_bins, mean_bins, std_bins, sem_bins, mean_abs_error, std_abs_error, sem_abs_error, mean_rel_error, std_rel_error, sem_rel_error] = binsDivision(realVec, estVec, numBins, minVec, maxVec, generationMinVal, generationMaxVal, param_name)
%binsDivision - Divide data to bins

% Screen needed range for parameters
general_good       = [];
good_idx           = [];

for i = 1:length(minVec)
    
    last_good_group    = good_idx;
    
    good_idx1          = find( realVec(i,:) >= minVec(i) );
    good_idx2          = find( realVec(i,:) <= maxVec(i) );
    good_idx           = intersect(good_idx1, good_idx2);
    
    general_good       = intersect(last_good_group, good_idx);
end

switch param_name
    case 'Flow'
        chosen_idx = 1;
        realData   = realVec(1,:);
        estData    = estVec(1,:);
        
    case 'Ktrans'
        chosen_idx = 2;
        realData   = realVec(2,:);
        estData    = estVec(2,:);
    case 'Vb'
        chosen_idx = 3;
        realData   = realVec(3,:);
        estData    = estVec(3,:);
    case 'Ve'
        chosen_idx = 4;
        realData   = realVec(4,:);
        estData    = estVec(4,:);
        
    case 'E'
        chosen_idx = 5;
        realData   = realVec(5,:);
        estData    = estVec(5,:);
        
    case {'BAT','Delay'}
        chosen_idx = 5; % Need to add delay min and max to .mat files
        realData   = realVec(6,:);
        estData    = estVec(6,:);
    otherwise
        error(['-E- Unknown parameter name: ' param_name ' !']);
end

real_data_filtered = realData(general_good);
est_data_filtered  = estData(general_good);

% Compute bins
binEdges = linspace( min(real_data_filtered), max(real_data_filtered), numBins + 1);
aj = binEdges(1:end-1);     %# bins lower edge
bj = binEdges(2:end);       %# bins upper edge
cj = ( aj + bj ) ./ 2;      %# bins center

% Assign values to bins
[~, binIdx] = histc(real_data_filtered, [binEdges(1:end-1) Inf]);

% Count number of values in each bin
nj = accumarray(binIdx', 1, [numBins 1], @sum);

% figure;
% bar(cj,nj,'hist')
% set(gca, 'XTick',binEdges, 'XLim',[binEdges(1) binEdges(end)])
% xlabel('Bins'), ylabel('Counts'), title('Histogram of Bins')

mean_bins  = zeros(numBins, 1);
std_bins   = zeros(numBins, 1);
sem_bins   = zeros(numBins, 1);

for i = 1 : numBins
    
    if sum(binIdx == i) > 0 % If any such exist
        relevant_est_data_idx = binIdx==i;    
        relevant_est_data     = est_data_filtered(relevant_est_data_idx);
        mean_bins(i)          = mean(relevant_est_data);
        std_bins(i)           = std(relevant_est_data);
        sem_bins(i)           = std(relevant_est_data) / sqrt(length(relevant_est_data));
    end
    
end

% Bound std values so mean+- std will not exceed range of parameter
upper_bound = mean_bins + abs(std_bins);
lower_bound = mean_bins - abs(std_bins);

% if isinf( minVec(chosen_idx) )
%     lowerLimit     = min(real_data_filtered);
% else
%     lowerLimit     = minVec(chosen_idx);
% end
% 
% if isinf( maxVec(chosen_idx) )
%     upLimit     = max(real_data_filtered);
% else
%     upLimit     = maxVec(chosen_idx);
% end

lowerLimit = generationMinVal(chosen_idx);
upLimit    = generationMaxVal(chosen_idx);
%upLimit    = min( max(real_data_filtered), generationMaxVal(chosen_idx));


upper_delta = max((upper_bound - upLimit)   , 0);
lower_delta = min((lower_bound - lowerLimit), 0);
std_bins    = std_bins - max(abs(upper_delta),abs(lower_delta));


% Assign the new real value vector, the middle of the bins
real_bins = cj';

abs_error = abs(   real_data_filtered - est_data_filtered);
rel_error = 100 * abs( ( real_data_filtered - est_data_filtered) ./ real_data_filtered ); % In percent

mean_abs_error = mean(abs_error);
std_abs_error  = std(abs_error);
sem_abs_error  = std(abs_error) / sqrt(length(abs_error));
mean_rel_error = mean(rel_error);
std_rel_error  = std(rel_error);
sem_rel_error  = std(rel_error) / sqrt(length(abs_error));

end

