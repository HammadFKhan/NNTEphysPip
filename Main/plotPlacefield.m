function Spikes = plotPlacefield(Spikes)

%Parse
figure('Name','Place Field Candidates');
trackLength = 200;
disp(['Track Length set to ' num2str(trackLength)])
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
            allplaceFields(trial,:) = Smooth(spikeRate,1);
            placeFieldperTrial(neuron,:,trial) = Smooth(spikeRate,3); % neuron x position x trial
        end
    end
    allplaceFieldsAvg(neuron,:) = mean(allplaceFields(1:8,:),1); %Change the numbers for trial
    position = 1:trackLength;
    s = ceil(sqrt(size(Spikes.PlaceFields.placeField,2)));
    subplot(s,s,neuron),spike_map(placeFieldMap(:,:,neuron),position);
    axis tight, axis off, colorbar off;
    title(['PF ' num2str(neuron)]);
end

for perTrial = 1:size(placeFieldperTrial,3)
    figure('Name','Place Field per Trial'),spike_map(placeFieldperTrial(:,:,perTrial),position);
end

figure('Name','All Place Fields across Trials')
spike_map(allplaceFieldsAvg,1:size(allplaceFieldsAvg,2));
Spikes.allplaceFields = allplaceFieldsAvg;
% Sort Place Fields
[sortedPlaceFields,normPlaceFields] = placefieldSort(allplaceFieldsAvg);
figure('Name','Sorted Place Fields'),spike_map(sortedPlaceFields,position);
figure('Name','Normalized Place Fields'),spike_map(normPlaceFields(:,1:100),position);
Spikes.PlaceFields.normPlaceFields = normPlaceFields;
end
