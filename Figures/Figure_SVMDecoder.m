%% Bulk processing and formating for svm decoder 

files = dir(fullfile('Y:\Hammad\Ephys\SeqProject\ForceField\','*.mat'));
for fileNum = 1:length(files)
    fName = fullfile(files(fileNum).folder,files(fileNum).name);
    disp(['Loading ' fName '...'])
    load(fName)
    if ~exist('fpath','var')
        [fpath,fname] = fileparts(ds_filename);
        error('No fpath detected!')
    end
    Spikes = leverPSTHSq(Spikes,IntanBehaviour);
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
    Spikes = rejectSpikes(Spikes,0.1,0.15,IntanBehaviour.parameters); % Reject spikes here for further analysis
    IntanBehaviour.reactionTime = 0;

    %Neural Trajectory Segementation using GPFA
    % Note that we concatenate trial conditions as to apply the same models for
    % statistical comparison (ie. hit vs miss, hit vs FA, opto vs no opto)
    Spikes = makeSpikeGPFA(Spikes);
    Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
    for n = 1:length(IntanBehaviour.hitTrace)%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
        Spikes.GPFA.HitMiss.dat(n).trialId = n;
    end
    Spikes.GPFA.MIHitFA.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.MIFA.dat];
    for n = length(IntanBehaviour.MIHitTrace)+1:length(Spikes.GPFA.MIHitFA.dat) %fix trials
        Spikes.GPFA.MIHitFA.dat(n).trialId = n;
    end
    Spikes.GPFA.HitEffort.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.effortperturb.dat];
    for n = 1:length(IntanBehaviour.MIHitTrace)%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
        Spikes.GPFA.HitEffort.dat(n).trialId = n;
    end
    addpath(genpath('C:\Users\khan332\Documents\GitHub\NeuralTraj'));
    addpath(genpath('mat_results'));
    if exist('mat_results','dir'),rmdir('mat_results','s'),end
    [Spikes.GPFA.resultHit,Spikes.GPFA.seqTrainHit] = gpfaAnalysis(Spikes.GPFA.hit.dat,1); %Run index
    [Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,2); %Run index
    [Spikes.GPFA.resultMIHit,Spikes.GPFA.seqTrainMIHit] = gpfaAnalysis(Spikes.GPFA.MIHit.dat,3); %Run index
    % [Spikes.GPFA.resultMIFA,Spikes.GPFA.seqTrainMIFA] = gpfaAnalysis(Spikes.GPFA.MIFA.dat,4); %Run index
    [Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,5); %Run index
    % [Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,6); %Run index
    [Spikes.GPFA.resultHitEffort,Spikes.GPFA.seqTrainHitEffort] = gpfaAnalysis(Spikes.GPFA.HitEffort.dat,7); %Run index
    close all

    % Now lets move all of this data into a new folder insider the
    % collected spikes directory for the rslds model to access which we can
    % call warpedSpks_sessions.
    targetDir = files(fileNum).folder;
    % Make a new directory folder if it does not exist
    newFolderName = 'svmSpks_sessions';
    fpath = fullfile(targetDir, newFolderName);
    % We create a new fpath so that rslds can reach it

    if ~exist(fpath, 'dir')
        mkdir(fpath);
        fprintf('Created new directory: %s\n', fpath);
    else
        fprintf('Directory already exists: %s\n', fpath);
    end
    sessionName = [fpath,'\',files(fileNum).name(1:end-4),'_warpedSpks.mat'];
    % save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
    save(sessionName,"Spikes","IntanBehaviour","ds_filename","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
    disp('Saved!')
    clear fpath fname sessionName warpedSpks Spikes IntanBehaviour ds_filename
end