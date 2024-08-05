%% Beta event intervals 
function BetaIEI(betaGroup,betaGroupRest)
figure,hold on
L23EventISI = [];
L23RunEventISI = [];
L23InitiateEventISI = [];
%L2/3
count = 1;
for ii = 1:30
    for i = 1:size(betaGroup(ii).electrode.betaBurst.NumDetectedBeta,1)
        if isempty(betaGroup(ii).electrode.betaBurst.detectedBeta{i})
            t = 0;t1 = 0;
        else
            t = betaGroup(ii).electrode.betaBurst.detectedBeta{i}(:,2);
            t1 = betaGroup(ii).electrode.betaBurst.window(i,1);
            win = t-t1;
%             plot(win,i*ones(1,length(t)),'.')
            % ERD (Event Rate Desynchonization)
            %         A = win(win<=1);
            %         B = win(win>1);
            %         ERD(i) = (length(A)-length(B))/(length(A)+length(B));
            %Beta Event Interval
            L23EventISI{count} = abs(diff(win));
            L23RunEventISI{count} = abs(diff(win(win>=1))); % Run Window
            L23InitiateEventISI{count} = abs(diff(win(win<1))); % Initiate Window
            if isempty(L23RunEventISI{count}) && isempty(L23InitiateEventISI{count})
                L23RunEventISI{count} = -.5;
                L23InitiateEventISI{count} = -.5;
            elseif isempty(L23RunEventISI{count}) 
                L23RunEventISI{count} = -.5;
            elseif isempty(L23InitiateEventISI{count})
                L23InitiateEventISI{count} = -.5;
            end
        end
        count = count+1;
    end
end

L23EventISI = vertcat(L23EventISI{:});
L23InitiateEventISI = vertcat(L23InitiateEventISI{:});
L23RunEventISI = vertcat(L23RunEventISI{:});

figure('Name','Layer 2/3'),histogram(L23EventISI,0:.1:2)
figure('Name','Initiate Layer 2/3'),yyaxis left,histogram(L23InitiateEventISI,0:.1:1),title('Initiate Layer2/3'),ylim([0 80])
yyaxis right,histogram(L23InitiateEventISI(L23InitiateEventISI==-0.5),-1:0.1:0),ylim([0 200])
figure('Name','Run Layer 2/3'),yyaxis left,histogram(L23RunEventISI,0:.1:1),title('Run Layer 2/3'),ylim([0 80])
yyaxis right,histogram(L23RunEventISI(L23RunEventISI==-0.5),-1:0.1:0),ylim([0 200])


% L5
count = 1;
figure,hold on
L5EventISI = [];
L5InitiateEventISI = [];
L5RunEventISI = [];
for ii = 30:64
    for i = 1:size(betaGroup(ii).electrode.betaBurst.NumDetectedBeta,1)
        if isempty(betaGroup(ii).electrode.betaBurst.detectedBeta{i})
            t = 0;t1 = 0;
        else
            t = betaGroup(ii).electrode.betaBurst.detectedBeta{i}(:,2);
            t1 = betaGroup(ii).electrode.betaBurst.window(i,1);
            win = t-t1;
%             plot(win,i*ones(1,length(t)),'.')
            % ERD (Event Rate Desynchonization)
            %         A = win(win<=1);
            %         B = win(win>1);
            %         ERD(i) = (length(A)-length(B))/(length(A)+length(B));
            L5EventISI{count} = abs(diff(win));
            L5RunEventISI{count} = abs(diff(win(win>=1))); % Run Window
            L5InitiateEventISI{count} = abs(diff(win(win<1))); % Initiate Window
            if isempty(L5RunEventISI{count}) && isempty(L5InitiateEventISI{count})
                L5RunEventISI{count} = -.5;
                L5InitiateEventISI{count} = -.5;
            elseif isempty(L5RunEventISI{count})
                L5RunEventISI{count} = -.5;
            elseif isempty(L5InitiateEventISI{count})
                L5InitiateEventISI{count} = -.5;
            end
        end
        count = count+1;
    end
end
L5EventISI = vertcat(L5EventISI{:});
L5InitiateEventISI = vertcat(L5InitiateEventISI{:});
L5RunEventISI = vertcat(L5RunEventISI{:});

figure('Name','Layer 5'),histogram(L5EventISI,0:.1:2)
figure('Name','Initiate Layer 5'),yyaxis left,histogram(L5InitiateEventISI,0:.1:1),title('Initiate Layer 5')%,ylim([0 45])
yyaxis right,histogram(L5InitiateEventISI(L5InitiateEventISI==-0.5),-1:0.1:0),ylim([0 200])
figure('Name','Run Layer 5'),yyaxis left,histogram(L5RunEventISI,0:.1:1),title('Run Layer 5')
yyaxis right,histogram(L5RunEventISI(L5RunEventISI==-0.5),-1:0.1:0),ylim([0 200])
% Rest State
count = 1;
figure,hold on
L23RestEventISI = [];
%L2/3
for ii = 30:64
    for i = 1:size(betaGroupRest(ii).electrode.betaBurst.NumDetectedBeta,1) % Keep the same number of trials
        if isempty(betaGroupRest(ii).electrode.betaBurst.detectedBeta{i})
            t = 0;t1 = 0;
        else
            t = betaGroupRest(ii).electrode.betaBurst.detectedBeta{i}(:,2);
            t1 = betaGroupRest(ii).electrode.betaBurst.window(i,1);
            win = t-t1;
            win = win(win<1); %% split in half to make the inwdow the same for run and intitate
            %             plot(win,i*ones(1,length(win)),'.')
            %Beta Event Interval
            L23RestEventISI{count} = abs(diff(win));
            if isempty(L23RestEventISI{count})
                L23RestEventISI{count} = -.5;
            end
        end
        count = count+1;
    end
end
L23RestEventISI = vertcat(L23RestEventISI{:});

figure('Name','Rest Layer 2/3'),yyaxis left,histogram(L23RestEventISI,0:.1:1),title('Rest Layer 2/3')
yyaxis right,histogram(L23RestEventISI(L23RestEventISI==-0.5),-1:0.1:0)%,ylim([0 350])
% L5
count = 1;
figure,hold on
L5RestEventISI = [];
for ii = 1:30
    for i = 1:size(betaGroupRest(ii).electrode.betaBurst.NumDetectedBeta,1)
        if isempty(betaGroupRest(ii).electrode.betaBurst.detectedBeta{i})
            t = 0;t1 = 0;
        else
            t = betaGroupRest(ii).electrode.betaBurst.detectedBeta{i}(:,2);
            t1 = betaGroupRest(ii).electrode.betaBurst.window(i,1);
            win = t-t1;
            win = win(win<1); %% split in half to make the inwdow the same for run and intitate
%             plot(win,i*ones(1,length(win)),'.')
            L5RestEventISI{count} = abs(diff(win));
            if isempty(L5RestEventISI{count})
                L5RestEventISI{count} = -.5;
            end
        end
        count = count+1;
    end
end
L5RestEventISI = vertcat(L5RestEventISI{:});
figure('Name','Rest Layer 5'),yyaxis left,histogram(L5RestEventISI,0:.1:1),title('Rest Layer 5'),ylim([0 3000])
yyaxis right,histogram(L5RestEventISI(L5RestEventISI==-0.5),-1:0.1:0),ylim([0 12000])
%%
figure,boxplot(L23InitiateEventISI(L23InitiateEventISI>0),'PlotStyle','compact'),ylim([0 1]),title('L23 Initiate')
figure,boxplot(L23RunEventISI(L23RunEventISI>0),'PlotStyle','compact'),ylim([0 1]),title('L23 Run')
figure,boxplot(L5InitiateEventISI(L5InitiateEventISI>0),'PlotStyle','compact'),ylim([0 1]),title('L5 Initiate')
figure,boxplot(L5RunEventISI(L5RunEventISI>0),'PlotStyle','compact'),ylim([0 1]),title('L5 Run')
%%
l23run = [length(L23RunEventISI(L23RunEventISI>0));length(L23RunEventISI(L23RunEventISI==-0.5))]/size(L23RunEventISI,1);
figure,pie(l23run),title('L23 Run')
legend('Burst','Single')

l23init = ([length(L23InitiateEventISI(L23InitiateEventISI>0));length(L23InitiateEventISI(L23InitiateEventISI==-0.5))])/size(L23InitiateEventISI,1);
figure,pie(l23init),title('L23 Initiate')
legend('Burst','Single')

l5run = [length(L5RunEventISI(L5RunEventISI>0));length(L5RunEventISI(L5RunEventISI==-0.5))]/size(L5RunEventISI,1);
figure,pie(l5run),title('L5 Run')
legend('Burst','Single')

l5init = [length(L5InitiateEventISI(L5InitiateEventISI>0));length(L5InitiateEventISI(L5InitiateEventISI==-0.5))]/size(L5InitiateEventISI,1);
figure,pie(l5init),title('L5 Initiate')
legend('Burst','Single')

l23rest = [length(L23RestEventISI(L23RestEventISI>0));length(L23RestEventISI(L23RestEventISI==-0.5))]/size(L23RestEventISI,1);
figure,pie(l23rest),title('L2/3 Rest')
legend('Burst','Single')

l5rest = [length(L5RestEventISI(L5RestEventISI>0));length(L5RestEventISI(L5RestEventISI==-0.5))]/size(L5RestEventISI,1);
figure,pie(l5rest),title('L5 Rest')
legend('Burst','Single')
%%
% %%
% figure,boxplot(L23EventISI,'PlotStyle','compact'),ylim([0 2]),title('L23')
% figure,boxplot(L5EventISI,'PlotStyle','compact'),ylim([0 2]),title('L5')
% figure,boxplot(L23RestEventISI,'PlotStyle','compact'),ylim([0 2]),title('L23 Rest')
% figure,boxplot(L5RestEventISI,'PlotStyle','compact'),ylim([0 2]),title('L5 Rest')
% %%
% figure,bar(sort(ERD)),ylim([-1 1])
% figure,bar(sort(ERDrest)),ylim([-1 1])
% figure, histogram(ERDrest,-1:0.5:1),hold on
% histogram(ERD,-1:0.5:1)
% % ERD probability
% ERDprop = length(ERD(abs(ERD)>=0.6))/length(ERD);
% ERDpropRest = length(ERDrest(abs(ERDrest)>=0.6))/length(ERDrest);
% figure,bar([ERDpropRest;ERDprop]),hold on
% % er = errorbar([ERDpropRest;ERDprop],std([ERDpropRest;ERDprop]'));    
% % er.Color = [0 0 0];                            
% % er.LineStyle = 'none';  