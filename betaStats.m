function stats = betaStats(bstats,LFPdepth)
[~,l23Idx] = find(LFPdepth<250);
[~,l4Idx] = find(LFPdepth>250 & LFPdepth<400);
[~,l5Idx] = find(LFPdepth>400);

LFPlayers = [repmat({'L23'},length(l23Idx),1);repmat({'L4'},length(l4Idx),1);repmat({'L5'},length(l5Idx),1)];
figure('Name','Beta Duration'),boxplot(bstats(:,1),LFPlayers,'plotstyle','compact');
figure('Name','Beta Amplitude'),boxplot(bstats(:,2),LFPlayers,'plotstyle','compact');
figure('Name','Beta Event Rate'),boxplot(bstats(:,3),LFPlayers,'plotstyle','compact');