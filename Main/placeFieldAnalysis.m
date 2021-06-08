function Spikes = placeFieldAnalysis(Spikes)
trackLength = 200;
trialStart = 1;
trialCutoff = 8;
disp(['Track Length set to ' num2str(trackLength) ' cm'])
disp(['Trial Cutoff set to ' num2str(trialCutoff) ' runs'])
pause(2)
disp('Starting analysis...')
for neuron = 1:size(Spikes.PlaceFields.placeField,2)
    for trial = 1:size(Spikes.PlaceFields.placeField,1)
        placeField = Spikes.PlaceFields.placeField{trial,neuron};
        spikeRate = zeros(trackLength,1);
        if isempty(placeField)
            placeField = NaN;
            placeFieldMap(trial,:,neuron) = zeros(1,trackLength);
            allplaceFields(trial,:) = zeros(1,trackLength);
            placeFieldperTrial(neuron,:,trial) = zeros(1,trackLength);
        else
            spikeRate(placeField(1):placeField(end)) = Spikes.VR(trial).spikeRate(placeField(1):placeField(end),neuron);
            placeFieldMap(trial,:,neuron) = Smooth(spikeRate,2); % trial x position x neuron
            allplaceFields(trial,:) = Smooth(spikeRate,2);
            placeFieldperTrial(neuron,:,trial) = Smooth(spikeRate,2); % neuron x position x trial
        end
    end
    allplaceFields(neuron,:) = mean(allplaceFields(trialStart:trialCutoff,:),1); %Change the numbers for trial
end

%% Sort Place Fields and normalize
FRmax = max(allplaceFields,[],2);
count = 1;
idx = zeros(size(allplaceFields,1),2);
for i = 1:size(allplaceFields,1)
    if isempty(find(allplaceFields(i,:)==FRmax(i,1), 1))
        idx(i,1) = i; % Index of neuron
        idx(i,2) = 0; % Position in space
        count = count+1;
    elseif FRmax(i,1) == 0
        count = count+1;
    else
        idx(count,1) = i;
        idx(count,2) = find(allplaceFields(i,:)==FRmax(i,1));
        count = count+1;
    end
end

% Now to sort the idx from least to greatest
[~,idx_sort] = sort(idx(:,2),1); % Ascending list of locations
sortedPlaceFields = allplaceFields(idx_sort,:);
padding = length(find(FRmax==0));
sortedPlaceFields = sortedPlaceFields(padding+1:end,:);

% Normalize Place Fields for visualization
normPlaceFields = zeros(size(sortedPlaceFields,1),size(sortedPlaceFields,2));
for ii = 1:size(sortedPlaceFields,1)
    maxV = max(sortedPlaceFields(ii,:));
    minV = min(sortedPlaceFields(ii,:));
    normPlaceFields(ii,:) = (sortedPlaceFields(ii,:)-minV)./(maxV-minV);
end
%% Population Correlation
% Calculate a �temporal scale factor�, which estimates an amount of time
% it took for population activity to become de-correlated (�life-time� of
% population activity). Firing rate histograms of all neurons within each
% trial type (see above) were binned into 100 msec bins. Thus, activity of
% the whole neuronal 4 population within each 100 msec time bin was
% described by a vector that had a length of the number of active neurons
% (population vector). Spearman rank correlation between all pairs of
% population vectors characterizing all time bins was calculated.
[rho,~] = corr(allplaceFields,'Type','Spearman');
rho(isnan(rho)) = 0;
%% Outputs
Spikes.PlaceFields.placeFieldMap = placeFieldMap; % PF candidates over trials
Spikes.PlaceFields.placeFieldperTrial = placeFieldperTrial; % Place fields of all neurons over trials (3 dimensions)
Spikes.PlaceFields.allplaceFields = allplaceFields; % Averaged place fields across trials 
Spikes.PlaceFields.sortedPlaceFields = sortedPlaceFields; % Sorted place fields based on FR
Spikes.PlaceFields.normPlaceFields = normPlaceFields; % Normalized by peak FR (visual purposes)
Spikes.PlaceFields.populationVector = rho;
disp('Done!')
