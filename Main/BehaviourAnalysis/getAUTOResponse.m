%% Plot Automatic chart
function [trialMask] = getAUTOResponse(Behaviour,expFlag)
% Behavior is only automatic sequence data;
% expFlag is if we pass through the IntanBehaviour flag
if ~exist('expFlag','var')
    expFlag = 0;
end

if expFlag
    sqPullMask = zeros(length(Behaviour.hitTrace),length(Behaviour.hitTrace(1).time));
    for n = 1:length(Behaviour.hitTrace)
        sqPullMask(n,Behaviour.hitTrace(n).pullCount) = 1;
    end
else
    [Behaviour.cleanedpullCounts, Behaviour.pullIndices] = cleanTimeoutSequences(Behaviour,0);
    sqPullMask = zeros(size(Behaviour.cleanedpullCounts));

    for n = 1:length(Behaviour.pullIndices)
        sqPullMask(n,Behaviour.pullIndices{n}) = 1;
    end
end
% Plot out interleaved flagging of the contextual cue data
% Return interleaved trials index of when each cue was initiated
trialMask = sqPullMask;
trialMask(trialMask==0) = NaN;
time = linspace(-Behaviour.parameters.windowBeforePull,Behaviour.parameters.windowAfterPull,size(trialMask,2));
figure
color = [46,49 149]/255;
for n = 1:size(trialMask,1)
    scatter(time,n*trialMask(n,:),'filled','MarkerFaceColor',color),hold on
end
xlim([-2.5, Behaviour.parameters.windowAfterPull])
ylim([0 16.5])
xline(0,'k','reward')

% pullIPI = cellfun(@(x) diff(x),Behaviour.pullIndices,'UniformOutput',false);
% pullIPI = vertcat(pullIPI{:});
% figure,
% histogram(pullIPI(:,1)/100,0:0.1:1.2,'Normalization','probability','EdgeColor','none'),hold on
% histogram(pullIPI(:,2)/100,0:0.1:1.2,'Normalization','probability','EdgeColor','none')
% xlabel('IPI (s)')
% ylabel('Probability')
% box off,set(gca,'tickdir','out')
% axis square