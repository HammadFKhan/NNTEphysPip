function Spikes = detectPlacefield(Spikes)

% Place fields were identified based on the firing rate of pyramidal cells.
% The mean firing rate as a function of position was computed for each cell
% in each trial condition. Regions on the track where the firing rate was
% above 20% of the peak were isolated. The length of these regions had to
% be longer than ?1/15th the length of the track and smaller than
% five-eighths the length of the track. The place cell also had to spike at
% least once while the subject was in the field on at least four-fifths of
% the trials.


for trial = 1:length(Spikes.VR)
    spikeRate = Spikes.VR(trial).spikeRate;
    for neuron = 1:size(spikeRate,2)
        mean_spikeRate = mean(spikeRate(:,neuron));
        rateThres = 1.2*mean_spikeRate; %Theshold set for rate
        thresholded = spikeRate(:,neuron)>rateThres;
        if sum(thresholded) < 2 %Need at least two values for start and stop
            disp(['Thresholding Failed for neuron ' num2str(neuron)]);
            Spikes.PlaceFields.placeField{trial,neuron} = [];
        else
            start = find(diff(thresholded)>0);
            stop = find(diff(thresholded)<0);
            if length(stop) == length(start)-1
                start = start(1:end-1);
            end
            % Exclude first place field if it is incomplete
            if length(stop)-1 == length(start)
                stop = stop(2:end);
            end
            % Correct special case when both first and last place fields are incomplete
            if start(1) > stop(1)
                stop(1) = [];
                start(end) = [];
            end
            firstPass = [start,stop];
            if isempty(firstPass)
                disp('Detection by thresholding failed')
                Spikes.PlaceFields.placeField{trial,neuron} = [];
            else
                disp(['After detection by thresholding: ' num2str(length(firstPass)) ' events.']);
                
                % Discard place fields that are too large
                placeField = [firstPass(:,1) firstPass(:,2)];
                duration = firstPass(:,2)-firstPass(:,1);
                max_placeField = .35*size(spikeRate(:,neuron),1);
                placeField(duration>max_placeField,:) = [];
                disp(['After max duration test: ' num2str(size(placeField,1)) ' events.']);
                
                % Discard place fields that are too small
                duration = firstPass(:,2)-firstPass(:,1);
                min_placeField = .03*size(spikeRate(:,neuron),1);
                placeField(duration<min_placeField,:) = [];
                disp(['After min duration test: ' num2str(size(placeField,1)) ' events.']);
                
                if isempty(placeField)
                    disp(['No place field detected for neuron ' num2str(neuron)]);
                    Spikes.PlaceFields.placeField{trial,neuron} = [];
                else
                    Spikes.PlaceFields.placeField{trial,neuron} = placeField;
                end
            end
        end
    end
end

end


