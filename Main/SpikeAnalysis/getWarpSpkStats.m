function warpedSpks = getWarpSpkStats(warpedSpks)
% I want to plot out the peak z scored value at the time of pull for each
% unit. To do this, we should calculate calculate the residuals of each
% unit as a funciton of the pulls. That will tell us which neurons is
% tuned to the response... I think
% Parameters
% Parameters
bin_size = 20; % example bin size - adjust to your data/time resolution
analysis_window = [-0.5, 0.5]; % analysis window in binned time indices (adjust as needed)
all_pulls = warpedSpks.warpSpikes;
% Get sizes from original data
[num_trials, num_timebins, num_neurons, num_pulls] = size(all_pulls);

% Calculate the number of bins after binning
n_time_binned = floor(num_timebins / bin_size);
binTime = linspace(warpedSpks.warpedTime(1,3),warpedSpks.warpedTime(end,3),n_time_binned);
% Preallocate binned data array
binned_all_pulls = zeros(n_time_binned, num_neurons, num_pulls);
disp('Applying bin...')
% Apply binning per neuron and pull
for p = 1:num_pulls
    for n = 1:num_neurons
        data = squeeze(all_pulls(:, :, n, p)); % trials × time
        binned_data = sum(getBin(data,bin_size))*(1000/bin_size);  % trials × binned_time
        binned_all_pulls(:, n, p) = binned_data;
    end
end

mean_baseline_neuron = squeeze(mean(binned_all_pulls(binTime<-2,:,:),[1,3]));
%%
% Preallocate mean response matrix: neurons × pulls
mean_responses = zeros(num_neurons, num_pulls);
wTime = warpedSpks.warpedTime;
binnedSpk = [];
% Calculate mean firing rate within analysis window per neuron/pull
% We want to normalize by the mean baseline
disp('Calculate mean response...')
for pull = 1:num_pulls
    win = find(wTime(:,pull)>=analysis_window(1) & wTime(:,pull)<=analysis_window(2));
    data = squeeze(all_pulls(:, win, :, pull)); % trials × bins × neurons
    % Here we calculate the binned z score response of the neurons
    for neuron = 1:num_neurons
        spkTemp = squeeze(data(:,:,neuron));
        binnedSpk(:,:,neuron) = getBin(spkTemp,bin_size);
    end
    avg_window = sum(binnedSpk,1);    % average over time bin dimension -> trials × 1 × neurons
    avg_window = squeeze(avg_window)*(1000/bin_size); % trials × neurons

    mean_responses(:, pull) = mean(avg_window, 1); % mean over trials, result is 1 × neurons
end
%%
%here we subtract the baseline response to calculate the residuals so we
%now know how responsive the neuron was to the stimulus
disp('Calculate stats...')
residuals = mean_responses-mean_baseline_neuron'; 
modulationIndex = (mean_responses-mean_baseline_neuron')./(mean_responses+mean_baseline_neuron');
selectivity_index = zeros(num_neurons, 1);
best_pull = zeros(num_neurons, 1);

for n = 1:num_neurons
    responses = abs(residuals(n, :));
    [R_best, idx_best] = max(responses);
    R_other = mean(responses(setdiff(1:num_pulls, idx_best)));
    responseTot(n,:) = responses;
    selectivity_index(n) = (R_best - R_other) / (R_best + R_other);
    best_pull(n) = idx_best;
end
%%% OUTPUT IT

warpedSpks.stats.residuals = residuals;
warpedSpks.stats.modulationIndex = modulationIndex;
warpedSpks.stats.selectivityIndex = selectivity_index;
warpedSpks.stats.bestPull = best_pull;

end