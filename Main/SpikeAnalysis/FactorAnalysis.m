load examgrades
[Loadings1,specVar1,T,stats] = factoran(grades,1);
stats.p
%%
[Loadings2,specVar2,T,stats] = factoran(grades,2,'rotate','none');
biplot(Loadings2, 'varlabels',num2str((1:5)'));
title('Unrotated Solution');
xlabel('Latent Factor 1'); ylabel('Latent Factor 2');
stats.dfe
%%
Sigma = cov(grades);
[LoadingsCov,specVarCov] = ...
        factoran(Sigma,2,'Xtype','cov','rotate','none');
%%
load carbig

X = [Acceleration Displacement Horsepower MPG Weight]; 
X = X(all(~isnan(X),2),:);

[Lambda,Psi,T,stats,F] = factoran(X,2,'Scores','regression');
inv(T'*T);   % Estimated correlation matrix of F, == eye(2)
Lambda*Lambda' + diag(Psi); % Estimated correlation matrix
Lambda*inv(T);              % Unrotate the loadings
F*T';                       % Unrotate the factor scores
%%
%Factor analysis was applied to spike counts in a 400-ms window, which
%either ended at stimulus onset (prestimulus) or began 100 ms after
%stimulus onset (stimulus). For the second V1 dataset (Fig. 6c), stimulus
%duration was only 400 ms, so we used a shorter (350 ms) window that began
%50 ms after stimulus onset. A data matrix, D, was compiled, with Dij being
%the spike count for the ith trial and the jth neuron. D included only
%neurons with rates ≥ 1 spike per s for both periods. Factor analysis was
%then applied to this matrix (neurons are variables, trials are
%observations). Factor analysis was performed separately on the D matrices
%compiled before and after stimulus onset. Factor analysis was performed
%separately for each condition. Subsequent results were averaged.

dat = PMDdata2(1:120);
Dpre = arrayfun(@(x) sum(x.spikes(:,1100:1600),2), dat, 'UniformOutput', false);
Dpost = arrayfun(@(x) sum(x.spikes(:,500:1000),2), dat, 'UniformOutput', false);
Dpre = horzcat(Dpre{:});
Dpost = horzcat(Dpost{:});
Dm1 = mean(Dpre,1);
Dm2 = mean(Dpost,1);
Dpre(:,Dm1<1 | Dm2<1) = [];
Dpost(:,Dm1<1 | Dm2<1) = [];

%% Check covariance
Dprecorr = corr(Dpre);
Dpostcorr = corr(Dpost);
figure
subplot(121),imagesc(Dprecorr);colormap(jet)
subplot(122),imagesc(Dpostcorr);colormap(jet)
%%
[Destpre,NetworkVarpre,PrivateNoisepre] = factorAnalysis(Dpre);
[Destpost,NetworkVarpost,PrivateNoisepost] = factorAnalysis(Dpost);
figure
subplot(221),imagesc(Destpre);colormap(jet);colorbar
subplot(222),imagesc(PrivateNoisepre);colormap(jet)
subplot(223),imagesc(Destpost);colormap(jet);colorbar
subplot(224),imagesc(PrivateNoisepost);colormap(jet)
%% Plot as a function of mean matched spikes
Dm1 = mean(Dpre,1);
Dm2 = mean(Dpost,1);

%%% find neurons that have similair spike rates
prebinSpike = discretize(Dm1,0:0.5:max([Dm1 Dm2]));
postbinSpike = discretize(Dm2,0:0.5:max([Dm1 Dm2]));
for n = 1:length(1:0.5:11)
    NoiseCorrelation(1,n) = mean(mean(Destpre(:,prebinSpike==n)));
    NoiseCorrelation(2,n) = mean(mean(Destpost(:,postbinSpike==n)));
end
%%% plot
figure
plot(NoiseCorrelation(1,:),'k--','LineWidth',2),hold on
plot(NoiseCorrelation(2,:),'k','LineWidth',2)
set(gca,'TickDir','out'),set(gca,'fontsize',12),box off
legend('Pre','Post')
%%
function [Dest,NetworkVar,PrivateNoise] = factorAnalysis(D)
[Lambda,Psi,T,stats,F] = factoran(D,2,'Scores','regression');
inv(T'*T)   % Estimated correlation matrix of F, == eye(2)
Dest = Lambda*Lambda' + diag(Psi); % Estimated correlation matrix which make up Network variance + private noise
NetworkVar = Lambda*Lambda';
PrivateNoise = diag(Psi);

figure
biplot(Lambda, 'varlabels',num2str((1:size(D,2))'));
title('Unrotated Solution');
xlabel('Latent Factor 1'); ylabel('Latent Factor 2');
end
