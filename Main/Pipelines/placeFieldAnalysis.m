function Spikes = placeFieldAnalysis(Spikes)
trackLength = 200;
trialStart = 17;
trialCutoff = 20;
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
            %Normalize
            placeFieldMap(trial,:,neuron) = Smooth(spikeRate,2); % trial x position x neuron
            maxV = max(placeFieldMap(trial,:,neuron));
            minV = min(placeFieldMap(trial,:,neuron));
            placeFieldMap(trial,:,neuron) = (placeFieldMap(trial,:,neuron)-minV)./(maxV-minV);
            allplaceFields(trial,:) = Smooth(spikeRate,3);
            placeFieldperTrial(neuron,:,trial) = Smooth(spikeRate,3); % neuron x position x trial
            numPlaceFields(trial) = length(nonzeros(sum(placeFieldperTrial(:,:,trial),2)));
        end
    end
    allplaceFields(neuron,:) = mean(allplaceFields(trialStart:trialCutoff,:),1); %Change the numbers for trial
end

%% Sort Place Fields, find spatial index, and normalize
% Place Fields per Trial
normplaceFieldperTrial{1,size(placeFieldperTrial,3)} = [];
for trial = 1:size(placeFieldperTrial,3)
    singlePFmap = placeFieldperTrial(:,:,trial);
    FRmax = max(singlePFmap,[],2);
    % Start sorting
    count = 1;
    idx = zeros(size(singlePFmap,1),2);
    for i = 1:size(singlePFmap,1)
        if isempty(find(singlePFmap(i,:)==FRmax(i,1), 1))
            idx(i,1) = i; % Index of neuron
            idx(i,2) = 0; % Position in space
            count = count+1;
        elseif FRmax(i,1) == 0
            count = count+1;
        else
            idx(count,1) = i;
            idx(count,2) = find(singlePFmap(i,:)==FRmax(i,1));
            count = count+1;
        end
    end
    
    % Now to sort the idx from least to greatest
    [~,idx_sort] = sort(idx(:,2),1); % Ascending list of locations
    sortedPlaceFields = singlePFmap(idx_sort,:);
    padding = length(find(FRmax==0));
    sortedPlaceFields = sortedPlaceFields(padding+1:end,:);
    
    %Find the information rate of all cells in this trial
    FRmean = mean(sortedPlaceFields,2); % mean firing rate for the cell across trials
    for j = 1:size(FRmean,1)
        lamda = sortedPlaceFields(j,:)/FRmean(j,1); %lamda = FR_i/mean_FR
        lamda(isnan(lamda)) = 0;
        for jj = 1:size(sortedPlaceFields,2)
            occupancy(:,jj) = jj./size(sortedPlaceFields,2);
        end
        information(j,:) = occupancy.*(lamda.*log2(lamda));
    end
    information(isnan(information))= 0;
    information = sum((information),2);
    informationRate{trial} = information;
    clear information
    
    % Normalize Place Fields
    normPlaceFields = zeros(size(sortedPlaceFields,1),size(sortedPlaceFields,2));
    for ii = 1:size(sortedPlaceFields,1) %Sort for each neuron
        maxV = max(sortedPlaceFields(ii,:));
        minV = min(sortedPlaceFields(ii,:));
        normPlaceFields(ii,:) = (sortedPlaceFields(ii,:)-minV)./(maxV-minV);
        %Place Field Width per neuron per trial
        findplaceFieldWidth = find(normPlaceFields(ii,:)>0.3);
        if isempty(findplaceFieldWidth)
            placeFieldWidth(ii,trial) = 0;
        else
            stop = diff(findplaceFieldWidth)==1;
            stop = find(stop == 0);
            if isempty(stop)
                placeFieldWidth(ii,trial) = abs(findplaceFieldWidth(end)-findplaceFieldWidth(1));
                centeredPF = Smooth(normPlaceFields(ii,findplaceFieldWidth(1):findplaceFieldWidth(end)),4)';
            else
                placeFieldWidth(ii,trial) = abs(findplaceFieldWidth(stop(1))-findplaceFieldWidth(1));
                centeredPF = Smooth(normPlaceFields(ii,findplaceFieldWidth(1):findplaceFieldWidth(stop(1))),4)';
            end
            width = placeFieldWidth(ii,trial)+1;
            centeredPF = padarray(centeredPF,[0 200-width],0,'post');
            centeredPlaceFields(ii,:,trial) = centeredPF;
        end
    end
    normplaceFieldperTrial{trial} = normPlaceFields;
    avg_centeredPlaceFields(trial,:) = mean(centeredPlaceFields(:,:,trial),1);
end

% Mean Place Fields
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
% Calculate a “spatial scale factor”, which estimates an amount of time
% it took for population activity to become de-correlated (‘life-time’ of
% population activity). Firing rate histograms of all neurons within each
% trial type (see above) were binned into 100 msec bins. Thus, activity of
% the whole neuronal 4 population within each 100 msec time bin was
% described by a vector that had a length of the number of active neurons
% (population vector). Spearman rank correlation between all pairs of
% population vectors characterizing all time bins was calculated.
[rho,~] = corr(allplaceFields(:,(1:100)),'Type','Spearman');
rho(isnan(rho)) = 0;
% Sets a flag for trunated analysis
if size(rho,1)<size(allplaceFields,2)
    disp(['Spatial scale factor set to ' num2str(size(rho,1)) ' bins instead of ' num2str(trackLength) ' bins'])
end
%% Remapping candidates
% Check FR for each place cell across trials
for trial = 1:size(normplaceFieldperTrial,2)
    for neuron = 1:size(normplaceFieldperTrial{1,trial},1)
    FR(neuron,trial) = mean(normplaceFieldperTrial{1,trial}(neuron,:),2);
    end
end

% Spike Data
%% Outputs
Spikes.PlaceFields.placeFieldMap = placeFieldMap; % PF candidates over trials
Spikes.PlaceFields.placeFieldperTrial = placeFieldperTrial; % Place fields of all neurons over trials (3 dimensions)
Spikes.PlaceFields.allplaceFields = allplaceFields; % Averaged place fields across trials
Spikes.PlaceFields.sortedPlaceFields = sortedPlaceFields; % Sorted place fields based on FR
Spikes.PlaceFields.normplaceFieldperTrial = normplaceFieldperTrial; % Sorted place fields per trial
Spikes.PlaceFields.normPlaceFields = normPlaceFields; % Normalized by peak FR (visual purposes)
Spikes.PlaceFields.populationVector = rho;
Spikes.PlaceFields.informationRate = informationRate;
Spikes.PlaceFields.placeFieldWidth = placeFieldWidth;
Spikes.PlaceFields.centeredPlaceFields = centeredPlaceFields;
Spikes.PlaceFields.avg_centeredPlaceFields = avg_centeredPlaceFields;
disp('Done!')

