%% Sparse CC analysis
% Script to generate CCA analysis of neural trajectory during task.

% Following the standard approach in CCA, we identified two sets of loading
% vectors, {wi} and {vi}, termed here as CCA modes, each of which was an
% activity mode within one of the two neural ensembles (that is, with N1
% and N2 elements, respectively). The index i ∈ {1, 2, 3, ..., minimum(N1,
% N2)} denoted the individual modes, which we determined such that the
% projections of the neural activity fluctuations, X and Y, onto wi and vi,
% were maximally correlated between the two trajectories, subject to the
% normalization constraint. Given this normalization condition, the
% quantity ) equals the correlation coefficient of the activity modes,and
% in the two different brain areas. After finding the first CCA mode (i
% =1), we identified successive modes in an iterative manner. Specifically,
% for all previously identified CCA modes we removed the CCA fluctuations
% from X and Y. We applied equation (11) to the residuals and thereby
% identified a set of orthonormal fluctuation modes with correlation
% coefficient values that progressively declined with the index, i. To
% identify the maxima specified by equation (11), we first randomly
% initialized the vectors wi and vi while constraining them to have unity
% length. We then found values of wi and vi that maximized the objective
% function in equation (11) by performing an alternating optimization

% Ebrahimi, S., Lecoq, J., Rumyantsev, O. et al. Emergent reliability in
% sensory cortical coding and inter-area communication. Nature 605, 713–721
% (2022).

% Load data 
if ~exist('M1Spikes','var')||~exist('M2Spikes','var')
    disp('Loading data...')
    load('Y:\Hammad\Ephys\LeverTask\DualShank\075356DualShank\Day3\Day3M2DualM1SingleRecording1_240727_180222\CCA_data')
end
% Make M1 and M2 PCA dimensions based on GPFA
[M1rh,M1rm,M1rmh,M1rmf] = trajNorm(M1Spikes,IntanBehaviour);
[M2rh,M2rm,M2rmh,M2rmf] = trajNorm(M2Spikes,IntanBehaviour);


%% Sparse CCA Analysis
% Here we take the high dimensional neural trajectory data and perform CCA
% analysis on it to see what correlations there are from the time varying
% signals. I chose the top 5 modes based on trajectories that occupy 15
% latent dimensions. 
% For statistical analysis we build seperate CCA models on subset of trials
% lets say we only use 80% of the data to check for validity.

nModes = 5;
iter = 2; %Number of training rounds
dataKeep = 0.8; % Percentage we keep for CCA model
timeLag = 10; % Control condition where we set the time lag for CCA control

CCA_timeLag10 = getCCA(M1rh,M1rm,M1rmh,M1rmf,M2rh,M2rm,M2rmh,M2rmf,iter,nModes,dataKeep,timeLag);


%%
figure,
subplot(121),imagesc(CCA.hit.rVec-mean(mean(CCA.hit.rVec))),hold on,colormap(jet)
subplot(122),imagesc(CCA.miss.rVec-mean(mean(CCA.hit.rVec))),colormap(jet),hold on
%% Plot average of Mode 1 to 3
f = figure;

CCAtype = CCA_timeLag2;

for n = 1:3
    dat = arrayfun(@(x) x.rVec(n,:),CCAtype.hit,'UniformOutput',false);
    dat = vertcat(dat{:});
    subplot(3,1,n),plot(smoothdata(mean(dat),'gaussian',25),'b'),hold on
    plot(smoothdata(mean(dat)+std(dat)/sqrt(iter),'gaussian',25),'b')
    plot(smoothdata(mean(dat)-std(dat)/sqrt(iter),'gaussian',25),'b')
    
    dat = arrayfun(@(x) x.rVec(n,:),CCAtype.miss,'UniformOutput',false);
    dat = vertcat(dat{:});
    plot(smoothdata(mean(dat),'gaussian',25),'r'),hold on
    plot(smoothdata(mean(dat)+std(dat)/sqrt(iter),'gaussian',25),'r')
    plot(smoothdata(mean(dat)-std(dat)/sqrt(iter),'gaussian',25),'r')
    box off, set(gca,'tickdir','out','fontsize',16),xlabel('Time'),ylabel('CC Coefficient'),axis square
end
f.Position = [681 159 560 800];

f = figure;

for n = 1:3
    dat = arrayfun(@(x) x.rVec(n,:),CCAtype.MIhit,'UniformOutput',false);
    dat = vertcat(dat{:});
    subplot(3,1,n),plot(smoothdata(mean(dat),'gaussian',25),'b'),hold on
    plot(smoothdata(mean(dat)+std(dat)/sqrt(iter),'gaussian',25),'b')
    plot(smoothdata(mean(dat)-std(dat)/sqrt(iter),'gaussian',25),'b')
    
    dat = arrayfun(@(x) x.rVec(n,:),CCAtype.MIFA,'UniformOutput',false);
    dat = vertcat(dat{:});
    plot(smoothdata(mean(dat),'gaussian',25),'r'),hold on
    plot(smoothdata(mean(dat)+std(dat)/sqrt(iter),'gaussian',25),'r')
    plot(smoothdata(mean(dat)-std(dat)/sqrt(iter),'gaussian',25),'r')
    box off, set(gca,'tickdir','out','fontsize',16),xlabel('Time'),ylabel('CC Coefficient'),axis square
end
f.Position = [681 159 560 800];

%% LOCAL FUNCTIONS
function r = meanTraj(X,trials,components)
r = X(:,:,trials);
r = permute(r,[1 3 2]);
end

function [rh,rm,rmh,rmf] = trajNorm(Spikes,Behaviour)
X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainHitMiss,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajHitMiss = reshape(X,size(X,1),Spikes.GPFA.seqTrainHitMiss(1).T,[]);

X = neuralTrajHitMiss;
hittrials = 1:length(Behaviour.cueHitTrace);
misstrials = length(Behaviour.cueHitTrace)+1:size(X,3);
rh = meanTraj(X,hittrials,6); %trajectory variable and predefined conditional trial indexes
rm = meanTraj(X,misstrials,6); %trajectory variable and predefined conditional trial indexes

X = arrayfun(@(x) vertcat(x.xorth),Spikes.GPFA.seqTrainMIHitFA,'UniformOutput',false);
X = horzcat(X{:});
neuralTrajMIHitFA = reshape(X,size(X,1),Spikes.GPFA.seqTrainMIHitFA(1).T,[]);

X = neuralTrajMIHitFA;

hittrials = 1:length(Behaviour.MIHitTrace);
FAtrials = length(Behaviour.MIHitTrace)+1:size(X,3);

rmh = meanTraj(X,hittrials,6);
rmf = meanTraj(X,FAtrials,6);

end

function CCA = getCCA(M1rh,M1rm,M1rmh,M1rmf,M2rh,M2rm,M2rmh,M2rmf,iter,nModes,dataKeep,timeLag)
% Initialize Structure
CCA = struct();
CCA.hit.wxMat = [];     CCA.miss.wxMat = [];       CCA.MIhit.wxMat = [];    CCA.MIFA.wxMat = [];
CCA.hit.wyMat = [];     CCA.miss.wyMat = [];       CCA.MIhit.wyMat = [];    CCA.MIFA.wyMat = [];
CCA.hit.rVec = [];      CCA.miss.rVec = [];        CCA.MIhit.rVec = [];    CCA.MIFA.rVec = [];

for nn = 1:iter
    datIdx = randperm(size(M1rh,2));
    datIdx = datIdx(1:floor(length(datIdx)*dataKeep));
    Xh = squeeze(M1rh(:,sort(datIdx),:));
    Yh = squeeze(M2rh(:,sort(datIdx),:));
    
    datIdx = randperm(size(M1rm,2));
    datIdx = datIdx(1:floor(length(datIdx)*dataKeep));
    Xm = squeeze(M1rm(:,sort(datIdx),:));
    Ym = squeeze(M2rm(:,sort(datIdx),:));
    
    datIdx = randperm(size(M1rmh,2));
    datIdx = datIdx(1:floor(length(datIdx)*dataKeep));
    Xmh = squeeze(M1rmh(:,sort(datIdx),:));
    Ymh = squeeze(M2rmh(:,sort(datIdx),:));
    
    datIdx = randperm(size(M1rmf,2));
    datIdx = datIdx(1:floor(length(datIdx)*dataKeep));
    Xmf = squeeze(M1rmf(:,sort(datIdx),:));
    Ymf = squeeze(M2rmf(:,sort(datIdx),:));
    
    for n = 1:size(M1rh,3)
        if timeLag %Checks if we want to do time lags
            X = squeeze(Xh(:,:,n));
            if (n+timeLag)<size(M1rh,3)
                Y = squeeze(Yh(:,:,n+timeLag)); %lead M2 by a certain amount
            else
                Y = squeeze(Yh(:,:,n)); %unless we reach the end of the timepoints
            end
        end
        [CCA.hit(nn).wxMat(:,:,n),CCA.hit(nn).wyMat(:,:,n),CCA.hit(nn).rVec(:,n)]=SparseCCA(X,Y,2,2,1,nModes);
        
        X = squeeze(Xm(:,:,n));
        Y = squeeze(Ym(:,:,n));
        [CCA.miss(nn).wxMat(:,:,n),CCA.miss(nn).wyMat(:,:,n),CCA.miss(nn).rVec(:,n)]=SparseCCA(X,Y,2,2,1,nModes);
        
        X = squeeze(Xmh(:,:,n));
        Y = squeeze(Ymh(:,:,n));
        [CCA.MIhit(nn).wxMat(:,:,n),CCA.MIhit(nn).wyMat(:,:,n),CCA.MIhit(nn).rVec(:,n)]=SparseCCA(X,Y,2,2,1,nModes);
        
        X = squeeze(Xmf(:,:,n));
        Y = squeeze(Ymf(:,:,n));
        [CCA.MIFA(nn).wxMat(:,:,n),CCA.MIFA(nn).wyMat(:,:,n),CCA.MIFA(nn).rVec(:,n)]=SparseCCA(X,Y,2,2,1,nModes);
        disp(['Timestep ' num2str(n) ' on iteration ' num2str(nn) '...'])
    end
end
end
