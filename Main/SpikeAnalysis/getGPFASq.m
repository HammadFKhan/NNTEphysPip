function [neuralDynamics] = getGPFASq(warpSpikes,IntanBehaviour)
% function to take in warped spikes and calculate gpfa neural trajectory
% data. From here we return the statistical structure of the neural
% trajectory properties from established pipelines were we call them in

%% Neural Trajectory Segementation using GPFA
% Note that we concatenate trial conditions as to apply the same models for
% statistical comparison (ie. hit vs miss, hit vs FA, opto vs no opto)

% Here we need to wrangle the data to import correctly. 

% warpSpikes is a 3 dimension array of trials x time x neurons from the
% warpedSpks structure that is curated during affine warp. 

tempSpikes = struct();

% Wrangle the warped data where the spks is a cell array of length neuron
% with each cell as trials x time
tempspks = cell(1,size(warpSpikes,3));
for n = 1:size(warpSpikes,3)
    tempspks{n} =  warpSpikes(:,:,n);
end
tempSpikes.PSTH.hit.spks = tempspks;

Spikes = makeSpikeGPFA(tempSpikes);
Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
for n = 1:length(IntanBehaviour.hitTrace)%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
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
% [Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,2); %Run index
% [Spikes.GPFA.resultMIHit,Spikes.GPFA.seqTrainMIHit] = gpfaAnalysis(Spikes.GPFA.MIHit.dat,3); %Run index
% [Spikes.GPFA.resultMIFA,Spikes.GPFA.seqTrainMIFA] = gpfaAnalysis(Spikes.GPFA.MIFA.dat,4); %Run index
% [Spikes.GPFA.resultHitMiss,Spikes.GPFA.seqTrainHitMiss] = gpfaAnalysis(Spikes.GPFA.HitMiss.dat,5); %Run index
% [Spikes.GPFA.resultMIHitFA,Spikes.GPFA.seqTrainMIHitFA] = gpfaAnalysis(Spikes.GPFA.MIHitFA.dat,6); %Run index
close all
%% Neural Trajectory Analysis
%IntanBehaviour.parameters = parameters;
%neuralTrajAnalysis(Spikes,Waves1,IntanBehaviour1);
[neuralDynamics,waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);
end