%% rslds prep only (no warping)
% Prep data for rslds without warping, useful when we want to include other
% trials types with varying pulls. For trials that are miss or FA with a
% non matching sq pull (for sample doesn't match 3 pulls), we fill the
% space with nan values to handle later down the road.
% For the state space analysis we don't care about how the data is arranged
% which is nice for computing the data. The issue comes with analyzing data
% down the road. Luckly we do a bulk of the analysis on the first pull,
% which makes a lot of sense. 

% We need something like this to match the rslds model which will be done
% on seperate dimensions.
% Holy smokes we can just index the trial types along the fourth dimensions
% since the time dim and neurons are the same. The only issue is the trials
% tho.... So we can instead transform the warpSpikes as a cell array with
% each cell entry being the trial type which we can affix with a label. 
files = dir(fullfile('Y:\Hammad\Ephys\SeqProject\ForceField\','*.mat'));
for fileNum = 2:length(files)
    fName = fullfile(files(fileNum).folder,files(fileNum).name);
    disp(['Loading ' fName '...'])
    load(fName)
    % Check if 'Spikes' exists
    if ~exist('Spikes', 'var')
        warning('Variable ''Spikes'' not found in %s. Skipping file.', fName);
        continue;
    end
    hitSpikes = zeros(size(Spikes.PSTH.hit.spks{1},1),size(Spikes.PSTH.hit.spks{1},2),length(Spikes.PSTH.hit.spks));
    for neuron = 1:length(Spikes.PSTH.hit.spks)
        hitSpikes(:,:,neuron) = Spikes.PSTH.hit.spks{neuron};
    end

    MIHitSpikes = zeros(size(Spikes.PSTH.MIHit.spks{1},1),size(Spikes.PSTH.MIHit.spks{1},2),length(Spikes.PSTH.MIHit.spks));
    for neuron = 1:length(Spikes.PSTH.MIHit.spks)
        MIHitSpikes(:,:,neuron) = Spikes.PSTH.MIHit.spks{neuron};
    end

    effortperturbSpikes = zeros(size(Spikes.PSTH.effortperturb.spks{1},1),size(Spikes.PSTH.effortperturb.spks{1},2),length(Spikes.PSTH.effortperturb.spks));
    for neuron = 1:length(Spikes.PSTH.effortperturb.spks)
        effortperturbSpikes(:,:,neuron) = Spikes.PSTH.effortperturb.spks{neuron};
    end

    [c,allTrials] = sort_hit_effort(IntanBehaviour); % Chronological built the outcome of spikes
    hiteffortperturbSpikes = cat(1,MIHitSpikes,effortperturbSpikes); % Combine trial types
    hiteffortperturbSpikes = hiteffortperturbSpikes(c,:,:);

    rsldsSpikes.taskSpikes{1} = hitSpikes;% <- organized as trialsxtimexneuron
    rsldsSpikes.taskSpikes{2} = MIHitSpikes;
    rsldsSpikes.taskSpikes{3} = effortperturbSpikes;
    rsldsSpikes.taskSpikes{4} = hiteffortperturbSpikes;


    trialTimetot(:,1) = -IntanBehaviour.parameters.windowBeforePull*IntanBehaviour.parameters.Fs:IntanBehaviour.parameters.windowAfterPull*IntanBehaviour.parameters.Fs;
    trialTimetot(:,2) = -IntanBehaviour.parameters.windowBeforeMI*IntanBehaviour.parameters.Fs:IntanBehaviour.parameters.windowAfterMI*IntanBehaviour.parameters.Fs;
    trialTimetot(:,3) = -IntanBehaviour.parameters.windowBeforeMI*IntanBehaviour.parameters.Fs:IntanBehaviour.parameters.windowAfterMI*IntanBehaviour.parameters.Fs;

    rsldsSpikes.taskTime = trialTimetot/1000; %ms to seconds
    rsldsSpikes.taskLabel = [{'reward-align hit'},{'mov-align hit'},{'mov-align effort'},{'mov-align hit-effort sorted'}];
    pullIndex = Spikes.PSTH.hit.pl;
    rsldsSpikes.hit.pull1 = pullIndex(:,1);
    rsldsSpikes.hit.pull2 = pullIndex(:,2);
    rsldsSpikes.hit.pull3 = pullIndex(:,3);
    rsldsSpikes.hit.leverTraces = horzcat(IntanBehaviour.hitTrace.trace);

    pullIndex = Spikes.PSTH.MIHit.pl;
    rsldsSpikes.MIHit.pull1 = pullIndex(:,1);
    rsldsSpikes.MIHit.pull2 = pullIndex(:,2);
    rsldsSpikes.MIHit.pull3 = pullIndex(:,3);
    rsldsSpikes.MIHit.leverTraces = horzcat(IntanBehaviour.MIHitTrace.trace);

    pullIndex = Spikes.PSTH.effortperturb.pl;
    rsldsSpikes.effortperturb.pull1 = pullIndex(:,1);
    rsldsSpikes.effortperturb.pull2 = pullIndex(:,2);
    rsldsSpikes.effortperturb.pull3 = pullIndex(:,3);

    pullIndex = [Spikes.PSTH.MIHit.pl;Spikes.PSTH.effortperturb.pl];
    rsldsSpikes.hiteffortperturb.pull1 = pullIndex(c,1);
    rsldsSpikes.hiteffortperturb.pull2 = pullIndex(c,2);
    rsldsSpikes.hiteffortperturb.pull3 = pullIndex(c,3);

    % Create labels so we know what trials are what exactly
    assert(length(allTrials)==size(pullIndex,1))
    % For brevity I assign 0 or 1 for noneffort trials (1 if no extra effort
    % was needed). We also made a label just in case we forget
    rsldsSpikes.hiteffortperturb.isnoeffort = allTrials(2,:);
    rsldsSpikes.hiteffortperturb.isnoeffortlabel = repmat({'hit'},size(pullIndex,1),1);
    rsldsSpikes.hiteffortperturb.isnoeffortlabel(allTrials(2,:)==0) = {'effort'};
    % Add the Intan Behavior traces to make life easy
    leverTrace = [horzcat(IntanBehaviour.MIHitTrace.trace),horzcat(IntanBehaviour.effortperturbTrace.trace)];
    rsldsSpikes.hiteffortperturb.leverTraces = leverTrace(:,c);

    % Now lets move all of this data into a new folder insider the
    % collected spikes directory for the rslds model to access which we can
    % call warpedSpks_sessions.
    targetDir = files(fileNum).folder;
    % Make a new directory folder if it does not exist
    newFolderName = 'rsldsSpks_sessions';
    fpath = fullfile(targetDir, newFolderName);
    % We create a new fpath so that rslds can reach it

    if ~exist(fpath, 'dir')
        mkdir(fpath);
        fprintf('Created new directory: %s\n', fpath);
    else
        fprintf('Directory already exists: %s\n', fpath);
    end

    sessionName = [fpath,'\',files(fileNum).name(1:end-4),'_rsldsSpks.mat'];
    save(sessionName,"rsldsSpikes","Spikes","IntanBehaviour","ds_filename","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
    disp('Saved!')
end
