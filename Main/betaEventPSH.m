function betaEventPSH(betaGroup,Spikes)
%% L2/3
% Fs set at 1024 for beta LFP
for i = 1:size(betaGroup(1).electrode.betaBurst.window,1)   % trials
    count = 1; 
    for ii = 1:30 % L2/3 beta events electrode
        betaEventTemp = cell2mat(betaGroup(ii).electrode.betaBurst.detectedBeta(i)); % beta event series based on electrode and trial
        for j = 1:size(betaEventTemp,1) % interate over all detected beta events across the single window
            betaEventTempadjusted = betaEventTemp(j,(1:3))-betaGroup(ii).electrode.betaBurst.window(i,1);
            spikeTemp = sort(vertcat(Spikes.spikes.L23Run(i).Cell{:})); % Combine L2/3 population
            % Find neurons that fire during beta event
            spikeInBetaTemp = spikeTemp((spikeTemp>=(betaEventTempadjusted(:,2)-.075) & spikeTemp<=betaEventTempadjusted(:,2)+.075)); 
            spikeTriggeredBeta(i).L23.spike{count} = spikeInBetaTemp; %Spike trigger for each beta event (good electrode x beta number)
            spikeTriggeredBeta(i).L23.betaEvent(count,:) = betaEventTempadjusted; % Beta Event timestamps and power
            spikeTriggeredBeta(i).L23.betaEventWin(count,:) = [betaEventTempadjusted(:,2)-.075, betaEventTempadjusted(:,2)+.075]; % Beta Event timestamps and power

            %             start = betaEventTemp(j,1)*1024;
            %             stop = betaEventTemp(j,3)*1024;
            peak = betaEventTemp(j,2)*1024;
            win = 0.075*1024;
            betaLine = betaGroup(ii).electrode.beta_band(peak-win:peak+win);
%             plot(spikeTemp(spikeInBetaTemp),i*ones(length(spikeInBetaTemp)),'.'),hold on
            plot(betaLine),hold on
            spikeTriggeredBeta(i).L23.betaLFP(count,:) = betaLine; %Beta event matched to spike
            count = count+1; % we use count to bypass indexing for j (ie. we dont care which beta event this is in the window)
        end
    end
end
%
% figure,subplot(2,1,1),hold on
% for i = 1:35
%     temp = horzcat(i*zeros(1,(spikeTriggeredBeta(1).L23.betaEvent(1,2)-0.075*1024),spikeTriggeredBeta(1).L23.betaLFP(1,:));
%     plot((1:length(temp))/1024,temp)
% end
% subplot(2,1,2), hold on
% for i = 1:35
%     plot(spikeTriggeredBeta(1).L23.spike{i},i*ones(length(spikeTriggeredBeta(1).L23.spike{i})),'.')
% end

%%
figure,
count = 1;
totalSpikes = [];
for ii = 1:size(betaGroup(1).electrode.betaBurst.window,1) 
    for i = 1:size(spikeTriggeredBeta(ii).L23.betaLFP,1)
        temp = spikeTriggeredBeta(ii).L23.betaLFP(i,:);
        subplot(2,1,1),plot((1:length(temp))/1024,temp),xlim([0 0.16]),hold on
        spikeInBetaTemp = spikeTriggeredBeta(ii).L23.spike{i}-spikeTriggeredBeta(ii).L23.betaEventWin(i,1);
        subplot(2,1,2),plot(spikeInBetaTemp ,count*ones(length(spikeInBetaTemp),1),'.'),xlim([0 0.16]), hold on
        count = count+1;
        totalSpikes = vertcat(totalSpikes,spikeInBetaTemp);
    end
end
temp1 = [];
for i = 1:size(betaGroup(1).electrode.betaBurst.window,1) 
    temp1 = vertcat(temp1,spikeTriggeredBeta(i).L23.betaLFP);
end
temp1 = mean(temp1,1);
subplot(2,1,1),plot((1:length(temp1))/1024,temp1,'k','LineWidth',3)

figure,histogram(totalSpikes,0:0.005:.16) % Spike Distrubution (pseudo spike rate)

%% Layer 5
for i = 1:size(betaGroup(1).electrode.betaBurst.window,1)  % trials
    count = 1; 
    for ii = 30:64 % L5 beta events electrode
        betaEventTemp = cell2mat(betaGroup(ii).electrode.betaBurst.detectedBeta(i)); % beta event series based on electrode and trial
        for j = 1:size(betaEventTemp,1) % interate over all detected beta events across the single window
            betaEventTempadjusted = betaEventTemp(j,(1:3))-betaGroup(ii).electrode.betaBurst.window(i,1);
            spikeTemp = sort(vertcat(Spikes.spikes.L5Run(i).Cell{:})); % Combine L2/3 population
            % Find neurons that fire during beta event
            spikeInBetaTemp = spikeTemp((spikeTemp>=(betaEventTempadjusted(:,2)-.075) & spikeTemp<=betaEventTempadjusted(:,2)+.075)); 
            spikeTriggeredBeta(i).L5.spike{count} = spikeInBetaTemp; %Spike trigger for each beta event (good electrode x beta number)
            spikeTriggeredBeta(i).L5.betaEvent(count,:) = betaEventTempadjusted; % Beta Event timestamps and power
            spikeTriggeredBeta(i).L5.betaEventWin(count,:) = [betaEventTempadjusted(:,2)-.075, betaEventTempadjusted(:,2)+.075]; % Beta Event timestamps and power

            %             start = betaEventTemp(j,1)*1024;
            %             stop = betaEventTemp(j,3)*1024;
            peak = betaEventTemp(j,2)*1024;
            win = 0.075*1024;
            betaLine = betaGroup(ii).electrode.beta_band(peak-win:peak+win);
%             plot(spikeTemp(spikeInBetaTemp),i*ones(length(spikeInBetaTemp)),'.'),hold on
            plot(betaLine),hold on
            spikeTriggeredBeta(i).L5.betaLFP(count,:) = betaLine; %Beta event matched to spike
            count = count+1; % we use count to bypass indexing for j (ie. we dont care which beta event this is in the window)
        end
    end
end
%%
figure,
count = 1;
totalSpikes = [];
for ii = 1:size(betaGroup(1).electrode.betaBurst.window,1) 
    for i = 1:size(spikeTriggeredBeta(ii).L5.betaLFP,1)
        temp = spikeTriggeredBeta(ii).L5.betaLFP(i,:);
        subplot(2,1,1),plot((1:length(temp))/1024,temp),xlim([0 0.16]),hold on
        spikeInBetaTemp = spikeTriggeredBeta(ii).L5.spike{i}-spikeTriggeredBeta(ii).L5.betaEventWin(i,1);
        subplot(2,1,2),plot(spikeInBetaTemp ,count*ones(length(spikeInBetaTemp),1),'.'),xlim([0 0.16]), hold on
        count = count+1;
        totalSpikes = vertcat(totalSpikes,spikeInBetaTemp);
    end
end


temp1 = [];
for i = 1:size(betaGroup(1).electrode.betaBurst.window,1) 
    temp1 = vertcat(temp1,spikeTriggeredBeta(i).L5.betaLFP);
end
temp1 = mean(temp,1);
subplot(2,1,1),plot((1:length(temp1))/1024,temp1,'k','LineWidth',3)
figure,histogram(totalSpikes,0:0.005:.16) % Spike Distrubution (pseudo spike rate)
