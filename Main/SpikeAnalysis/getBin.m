function binned_data = getBin(data,bin_size)
% Bin spike data
[n_trials, n_time] = size(data);

% Trim time dimension to multiple of bin_size
n_time_trim = floor(n_time / bin_size) * bin_size;
data_trim = data(:, 1:n_time_trim);

% Reshape to [n_trials, bin_size, n_time_trim/bin_size]
data_reshaped = reshape(data_trim', bin_size, [], n_trials); 
% Note the transpose is so time is first dimension for reshaping

% Sum or average within bins (along first dimension)
binned_data = squeeze(mean(data_reshaped, 1))';  
% Output size: [n_trials, n_time_trim/bin_size]
end