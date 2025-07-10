function [r,zscore_s,stability,matrix] = getMeanTraj(X,trials,components)
matrix = X(1:components,:,trials);
r = squeeze(mean(matrix,3));
% Subtract the mean from each element
centered_matrix = matrix - r;

% Calculate the Euclidean distance for each row
s = squeeze(sqrt(sum(centered_matrix.^2, 1)));

% baseline_s = s(1:70, :, :); % Pre-cue period
%zscore_s = (s - mean(baseline_s(:))) ./ std(baseline_s(:));
zscore_s = (s - min(s(:))) ./ (max(s(:)) - min(s(:)));

stability = 1 ./ (s + eps);
% Alternative: Z-score stability relative to baseline period
baseline_stability = stability(1:70, :, :); % Pre-cue period
stability = (stability - mean(baseline_stability(:))) ./ std(baseline_stability(:));
end