function stats = betaStats(bstats,LFPdepth)
[~,l23Idx] = find(LFPdepth<400);
% [~,l4Idx] = find(LFPdepth>250 & LFPdepth<400);
[~,l5Idx] = find(LFPdepth>400);
% bstats = flip(bstats);
LFPlayers = [repmat({'L23'},length(l23Idx),1);repmat({'L5'},length(l5Idx),1)];
% figure('Name','Beta Duration'),boxplot(bstats(:,1),LFPlayers,'plotstyle','compact');title('Beta Duration')
figure('Name','Beta AmplitudeNorm'),boxplot(cell2mat(bstats(:,2)),LFPlayers,'plotstyle','traditional');title('Beta Amplitude Norm')
figure('Name','Beta Amplitude'),boxplot(cell2mat(bstats(:,4)),LFPlayers,'plotstyle','traditional');title('Beta Amplitude')
% figure('Name','Beta Event Rate'),boxplot(cell2mat(bstats(:,7)),LFPlayers,'plotstyle','traditional');title('Beta Event Rate')
% figure('Name','Amplitude Fano vs CV'),plot(bstats(:,10),sqrt(bstats(:,11)),'.'),axis([0 1.5 0 1.5]);xlabel('FanoFactor'),ylabel('CV'),title('Amplitude FanoFactor')
