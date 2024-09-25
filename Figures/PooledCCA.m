%% Pool CCA responses across sessions
clear
files = dir(fullfile('D:\M1M2DualShank\CCA\','*.mat'));
for fileNum = 1:length(files)
    fName = fullfile(files(fileNum).folder,files(fileNum).name);
    disp(['Loading ' fName '...'])
    load(fName)
    % Extract PA structure/values
    CCAall(fileNum).CCA = CCA;
end
%%
hitCCA = []; missCCA = []; MIhitCCA = []; MIFACCA = [];
for n = 1:length(CCAall)
    dat = arrayfun(@(x) x.rVec(1,:),CCAall(n).CCA.hit,'UniformOutput',false);
    dat = vertcat(dat{:});
    hitCCA = vertcat(hitCCA,abs(mean(dat(:,75:end),2)-mean(dat(:,1:74),2)));
    
    dat = arrayfun(@(x) x.rVec(1,:),CCA.miss,'UniformOutput',false);
    dat = vertcat(dat{:});
    missCCA = vertcat(missCCA,abs(mean(dat(:,75:end),2)-mean(dat(:,1:74),2)));
    
    dat = arrayfun(@(x) x.rVec(1,:),CCA.MIhit,'UniformOutput',false);
    dat = vertcat(dat{:});
    MIhitCCA = vertcat(MIhitCCA,abs(mean(dat(:,75:end),2)-mean(dat(:,1:74),2)));
    
    dat = arrayfun(@(x) x.rVec(1,:),CCA.MIFA,'UniformOutput',false);
    dat = vertcat(dat{:});
    MIFACCA = vertcat(MIFACCA,abs(nanmean(dat(:,75:end),2)-nanmean(dat(:,1:74),2)));
end


figure,
subplot(121),customBoxplot([hitCCA, missCCA])
box off, set(gca,'tickdir','out','fontsize',16),ylabel('CC Coefficient')
subplot(122),customBoxplot([MIhitCCA, MIFACCA])
box off, set(gca,'tickdir','out','fontsize',16),ylabel('CC Coefficient')
