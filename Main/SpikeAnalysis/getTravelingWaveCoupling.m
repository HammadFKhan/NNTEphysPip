%% Spikes analysis
[fpath,name,exts] = fileparts(ds_filename);
data = matfile(ds_filename);
path = [fpath,'/kilosort3/'];
mergename = 'merged';
Kilosort3AutoMergeTester
path = [fpath,'/kilosort3/' mergename];
%%% Spike preprocessing (includes merging (optional) and channel info
%%% return)
% Read in kilosort data for matlab analysis
SpikeClusters = readNPY(fullfile(path, 'spike_clusters.npy'));
SpikeSamples = readNPY(fullfile(path, 'spike_times.npy'));
SpikeChannel = readNPY(fullfile(path,'channel_positions.npy'));
Spikes.SpikeClusters = SpikeClusters; 
Spikes.SpikeSamples = SpikeSamples;
Spikes = clusterSort(Spikes);
Spikes = ISI(Spikes,0.01,data.Fs,0); %Spikes, Interval, Fs
% Calculate Depth profile
%load chanMap64F2
load chanMap64Sharp
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms, max_site] =...
    spikeTemplatePosition(data.fpath,ycoords,[]); % 'invert'
for i = 1:length(tempAmps)
    Spikes.Clusters(i).spikeDepth = templateDepths(i);
    Spikes.Clusters(i).channelDepth = max_site(i);
    Spikes.Clusters(i).spikeAmplitude = tempAmps(i);
    Spikes.Clusters(i).waveforms = waveforms(i,:);
    Spikes.Clusters(i).spikeDuration = templateDuration(i)/data.Fs*1000;
end
%%% delete empty spikes
temp = arrayfun(@(x) isempty(x.cluster), Spikes.Clusters);
Spikes.Clusters(temp) = []; 
%%% Calculate trial PSTH for lever
Spikes = leverPSTH(Spikes,IntanBehaviour);
%%% save spike output data to load into gui
savepath = fullfile(path,['spks4sorting','.mat']);
path = [fpath,'/kilosort3/' mergename];
save(savepath,'Spikes','-v7.3')
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

%%% Neural Trajectory Segementation using GPFA
% Note that we concatenate trial conditions as to apply the same models for
% statistical comparison (ie. hit vs miss, hit vs FA, opto vs no opto)
Spikes = makeSpikeGPFA(Spikes);
Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
for n = 1:IntanBehaviour.nCueHit%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
    Spikes.GPFA.HitMiss.dat(n).trialId = n;
end
Spikes.GPFA.MIHitFA.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.MIFA.dat];
for n = length(IntanBehaviour.MIHitTrace)+1:length(Spikes.GPFA.MIHitFA.dat) %fix trials
    Spikes.GPFA.MIHitFA.dat(n).trialId = n;
end
%%%
addpath(genpath('C:\Users\khan332\Documents\GitHub\NeuralTraj'));
addpath(genpath('mat_results'));
if exist('mat_results','dir'),rmdir('mat_results','s'),end
[Spikes.GPFA.resultHit,Spikes.GPFA.seqTrainHit] = gpfaAnalysis(Spikes.GPFA.hit.dat,1); %Run index
[Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,2); %Run index
[Spikes.GPFA.resultMIHit,Spikes.GPFA.seqTrainMIHit] = gpfaAnalysis(Spikes.GPFA.MIHit.dat,3); %Run index
[Spikes.GPFA.resultMIFA,Spikes.GPFA.seqTrainMIFA] = gpfaAnalysis(Spikes.GPFA.MIFA.dat,4); %Run index
[Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,5); %Run index
[Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,6); %Run index
close all


%%% Neural Trajectory Analysis
%IntanBehaviour.parameters = parameters;
%neuralTrajAnalysis(Spikes,Waves1,IntanBehaviour1);
[M1neuralDynamics,M1waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
[spath,~,~] = fileparts(fpath);
[~,sname,~] = fileparts(spath);
dse = 'Y:\Hammad\Ephys\LeverTask\Data_for_Figures\TrajectoryWaveCoupling';
sessionName = [dse,'\',sname(1:end-10),'_SpikeWaveCoupling.mat'];
try
    rmLFPIntanBehaviour
catch
end
save(sessionName,"M1neuralDynamics","Spikes","M1waveDynamics","IntanBehaviour","fpath","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",
disp('Saved!')