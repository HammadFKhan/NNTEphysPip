function output = genEigenDecomp(rawData,Rdata,Sdata)
% Generalized eigenvalue decomposition
% Provide raw data for weighted transform (recommended raw data for ease of
% analysis)
% Rdata is reference data for covaraince map and Sdata is reference sample
% of interest
showplot = 1;
ntrials = size(Rdata,3);
nbchans = size(Rdata,1);
[allCovS,allCovR] = deal(zeros(ntrials,nbchans,nbchans)); % Initiate matrices
% Build SxR channels pre- post- stimulus
for triali = 1:ntrials
    seg = Sdata(:,:,triali); % extract one trial
    seg = seg-mean(seg,2); % temporally mean center
    allCovS(triali,:,:) = seg*seg'/size(seg,2); % Compute covariance with weight factor N
    % Do the same for reference matrix
    seg = Rdata(:,:,triali); % extract one trial
    seg = seg-mean(seg,2); % temporally mean center
    allCovR(triali,:,:) = seg*seg'/size(seg,2); % Compute covariance with weight factor N
end
% Clean up covariances
% Clean R
meanR = squeeze(mean(allCovR)); % average across trials
dists = zeros(ntrials,1);
for segn = 1:size(allCovR,1)
    r = allCovR(segn,:,:);
    % Euclidean distance
    dists(segn) = sqrt(sum((r(:)-meanR(:)).^2)); %Summed difference
end
% Clean trial outliers with RHO value >0.001 (2.98 z)
covR = squeeze(mean(allCovR(zscore(dists)<3,:,:),1));

% Now we do the same with S
meanS = squeeze(mean(allCovS)); % average across trials
dists = zeros(ntrials,1);
for segn = 1:size(allCovS,1)
    r = allCovS(segn,:,:);
    % Euclidean distance
    dists(segn) = sqrt(sum((r(:)-meanS(:)).^2)); %Summed difference
end
% Clean channel outliers with RHO value >0.001 (2.98 z)
covS = squeeze(mean(allCovS(zscore(dists)<3,:,:),1));

% GED calculation
[evecs,evals] = eig(covS,covR);
[evals,sidx]  = sort(diag(evals),'descend');
evecs = evecs(:,sidx);
data2D = reshape(rawData,nbchans,[]); % Raw data is reshaped into just electrode x time
% Calculate components
compts = evecs(:,1:3)'*data2D; % "component electrode" x time series
compts = reshape(compts,[],size(rawData,2)); % change back into trials
%%% power spectrum
comppower = abs(fft(compts,[],1)).^2;
comppowerAve = squeeze(mean(comppower,2));

%%% component map
compmap = evecs(:,1:3)' * covS;
% flip map sign
[~,se] = max(abs( compmap ));
compmap = compmap .* sign(compmap(se));

%%% Outputs
output.covS = covS;
output.covR = covR;
output.evecs = evecs;
output.evals = evals;
output.compts = compts;
output.comppowerAve = comppowerAve;
output.compmap = compmap;

%%% plots data
if showplot
   figure,subplot(131),imagesc(covR),colormap(hot)
    axis square
    title('Covariance R')
    xlabel('Electrode'), ylabel('Electrode')
    v = caxis;
    subplot(132),imagesc(covS),colormap(hot)
    axis square
    title('Covariance S')
    xlabel('Electrode'), ylabel('Electrode')
    caxis(v)
    subplot(133),imagesc(covS-covR),colormap(hot)
    axis square
    title('S-R')
    xlabel('Electrode'), ylabel('Electrode')
    
    figure,plot(evals,'ks-','markersize',10,'markerfacecolor','r')
    axis square, box off,set(gca,'tickdir','out')
    set(gca,'xlim',[0 30])
    title('GED components plot')
    xlabel('Component number'), ylabel('Power ratio (\lambda)')
    figure,imagesc(evecs(1:3,:)'),colormap(hot)
    title('Component Electrode Map')
    xlabel('Component number'), ylabel('Electrode')
    figure,plot(smoothdata(abs(compmap(1:3,:))'),'LineWidth',2)
    cmap = gray(4);
    set(gca(),'ColorOrder',cmap)
    title('Components Depth')
    xlabel('Depth'), ylabel('Component Power')
    box off,set(gca,'TickDir','out'),set(gca,'FontSize',16)
    legend
    
    figure,plot(squeeze(mean(compts,3)'));
    set(gca(),'ColorOrder',cmap)
    title('Components Time-series')
    xlabel('Samples'), ylabel('Component Amplitude')
    box off,set(gca,'TickDir','out'),set(gca,'FontSize',16)
    legend
end
