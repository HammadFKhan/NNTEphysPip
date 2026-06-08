%% Batch processing of incomplete sequences across animals
% Grab behavioural incomplete sequences
% Pass the windows for spike sorting
% Save as new spike array
clear
clc
files = dir(fullfile('Y:\Hammad\Ephys\SeqProject\SqOnly\DLS','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    IntanBehaviour = getIncompleteSequences(IntanBehaviour);
    %%% Spikes analysis
    % Calculate Depth profile
    %load chanMap64F2
    load chanMap64Sharp
    %%% Calculate trial PSTH for lever
    Spikes = leverPSTHSq(Spikes,IntanBehaviour);
    %%% save spike output data to load into gui
    %ManualSpikeCurateGUI
    %%% Basic spike analysis
    % z-score spike rates
    if exist('parameters','var')
        IntanBehaviour.parameters = parameters;
    end
    if exist('goodSpkComponents','var')
        Spikes.goodSpkComponents = unique(goodSpkComponents);
    else
        Spikes.goodSpkComponents = 1:length(Spikes.Clusters);
    end
    Spikes = rejectSpikes(Spikes,0.25,0.25,IntanBehaviour.parameters); % Reject spikes here for further analysis
    [Spikes] = sortSpkLever(Spikes,IntanBehaviour);

    % save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
    sessionName = fullfile(fullfile(files(fileNum).folder,files(fileNum).name));
    disp(['Saving ' num2str(sessionName) '...'])
    save(sessionName,"Spikes","IntanBehaviour","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
    disp('Saved!')
end