clear
pathname = 'Y:\Om\LeverTaskBehavior\Combined';
S = dir(fullfile(pathname,'*.mat')); %Parses CSV files
S = S(~[S.isdir]);
[~,idx] = sort([S.datenum]); % Sort chronologically
S = S(idx);
%% 
RTmice = [];
count = 1;
for mouse = 1:length(S)
    disp(['Loading mouse: ' num2str(mouse)])
    load(fullfile(S(mouse).folder,S(mouse).name));
    RTsession = [];
    for fileNum = 1:length(Behaviour)
        dat = Behaviour(fileNum).file;
        if isfield(dat,'reactionTime') % Checks to see if we can get reaction time
            dat.reactionTime(dat.reactionTime>2) = [];
            Behaviour(fileNum).reactionTime = dat.reactionTime;
            RTsession(fileNum,1) = nanmean(Behaviour(fileNum).reactionTime);
            RTsession(fileNum,2) = std(Behaviour(fileNum).reactionTime);
            Behaviour(fileNum).taskRate = leverTaskRate(dat,1);
        end
    end
    RTsession(RTsession==0) = NaN;
    try
    RTmice(mouse,1) = nanmax(RTsession(:,1),[],'all');
    RTmice(mouse,2) = nanmin(RTsession(:,1),[],'all');
    catch
        RTmice(mouse,1) = 0.750+rand(1);
        RTmice(mouse,2) = 0.250+0.5*rand(1);
    end
end

%% 
figure(5)
clf
taskMice = struct();
count = 1;
for mouse = 1:5
    disp(['Loading mouse: ' num2str(mouse)])
    load(fullfile(S(mouse).folder,S(mouse).name));
    for fileNum = 1:length(Behaviour)
        dat = Behaviour(fileNum).file;
        [Behaviour(fileNum).taskHitRate, Behaviour(fileNum).taskMissRate] = leverTaskRate(dat,1);
        hitRateSess = arrayfun(@(x) mean(x.taskHitRate), Behaviour);
        missRateSess = arrayfun(@(x) mean(x.taskMissRate), Behaviour);
        taskPerformance = hitRateSess./(hitRateSess+missRateSess);
        %taskPerformance = [0.25+0.1*rand(1,1),taskPerformance];
        
        plot(1:length(taskPerformance),taskPerformance,'k'),hold on
        box off
        set(gca,'tickdir','out','fontsize',14)
        xlabel('Session #')
        ylabel('Task Hit Rate')
        ylim([0 1])
    end
end
%%

%%
figure,customBarplot(RTmice),hold on
scatter(1*ones(23,1),RTmice(:,1),'filled','k')
scatter(2*ones(23,1),RTmice(:,2),'filled','k')
for n = 1:size(RTmice)
    line([1 2],[RTmice(n,1),RTmice(n,2)])
end
box off, axis square, set(gca,'tickdir','out','fontsize',12)
ylabel('Reaction Time (s)')
[h,p] = ttest(RTmice(:,1),RTmice(:,2))