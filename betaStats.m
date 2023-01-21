function stats = betaStats(bstats,LFPdepth)
[~,l23Idx] = find(LFPdepth<400);
% [~,l4Idx] = find(LFPdepth>250 & LFPdepth<400);
[~,l5Idx] = find(LFPdepth>400);
% bstats = flip(bstats);
LFPlayers = [repmat({'L23'},length(l23Idx),1);repmat({'L5'},length(l5Idx),1)];
% figure('Name','Beta Duration'),boxplot(bstats(:,1),LFPlayers,'plotstyle','compact');title('Beta Duration')
figure('Name','Beta AmplitudeNorm'),boxplot(cell2mat(bstats(:,2)),LFPlayers,'plotstyle','traditional'),...
    title('Beta Amplitude Norm'),box off, set(gca,'TickDir','out');
figure('Name','Beta Amplitude'),boxplot(cell2mat(bstats(:,4)),LFPlayers,'plotstyle','traditional'),...
    title('Beta Amplitude'),box off, set(gca,'TickDir','out');
L23betaEvent = cell2mat(bstats(l23Idx,8));
L5betaEvent = cell2mat(bstats(l5Idx,8));
betaEventLayers = [repmat({'L23'},length(L23betaEvent),1);repmat({'L5'},length(L5betaEvent),1)];
figure('Name','Beta AmplitudeSingles'),boxplot(cell2mat(bstats(:,8)),betaEventLayers,'plotstyle','traditional'),...
    title('Beta Amplitude'),box off, set(gca,'TickDir','out');


sElectrode = 10; %superficial electrode for beta event rate
dElectrode = 60; %deep electrode for beta even rate
L23 = repmat({'Superficial'},length(cell2mat(bstats(sElectrode,7))),1);
L5 =  repmat({'Deep'},length(cell2mat(bstats(dElectrode,7))),1);

figure('Name','Beta Event Rate'),boxplot(cell2mat(bstats([sElectrode dElectrode],7)),[L23;L5],'plotstyle','traditional');title('Beta Event Rate')

L23_1 = repmat({'Superficial'},length(cell2mat(bstats(sElectrode,1))),1);
L5_1 =  repmat({'Deep'},length(cell2mat(bstats(dElectrode,1))),1);

figure('Name','Beta Event Duration'),boxplot(cell2mat(bstats([sElectrode dElectrode],1)),[L23_1;L5_1],'plotstyle','traditional');title('Beta Event Duration')
% figure('Name','Amplitude Fano vs CV'),plot(bstats(:,10),sqrt(bstats(:,11)),'.'),axis([0 1.5 0 1.5]);xlabel('FanoFactor'),ylabel('CV'),title('Amplitude FanoFactor')

stats.BetaAmplitudeNorm = cell2mat(bstats(:,2));
stats.BetaEventRate = cell2mat(bstats([sElectrode dElectrode],7));
stats.BetaEventRateLabel = [L23;L5];
stats.BetaEventDuration = cell2mat(bstats([sElectrode dElectrode],1));
stats.BetaEventDurationLabel = [L23_1;L5_1];