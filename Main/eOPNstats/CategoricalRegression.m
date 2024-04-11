%% Categorical regression of eOPN M2
load('Y:\Hammad\Ephys\LeverTask\eOPNM2\23949eOPNRbp4\Day4\baselineRTSpeed.mat')
Condition = repmat("baseline",length(PQ),1); % Make the categorical var
IndNeuralSpeed = PQ; % Investigate the effect of NT Speed on RT
DepReactionT  = rawreactionTime; % Output (RT)
load('Y:\Hammad\Ephys\LeverTask\eOPNM2\23949eOPNRbp4\Day4\eOPNRTSpeed.mat')
Condition = [Condition;repmat("eOPN",length(PQ),1)]; % Make the categorical var
IndNeuralSpeed = [IndNeuralSpeed;PQ]; % Investigate the effect of NT Speed on RT
DepReactionT  = [DepReactionT;rawreactionTime]; % Output (RT)
experiment = table(DepReactionT,IndNeuralSpeed,Condition); % Make a table of how reactiontime depends on neural speed
experiment.Condition = categorical(experiment.Condition);

fit = fitlm(experiment,'DepReactionT~IndNeuralSpeed*Condition')

w = linspace(min(IndNeuralSpeed),max(IndNeuralSpeed));
figure()
gscatter(IndNeuralSpeed,DepReactionT,Condition,'bgr','x.o')
line(w,feval(fit,w,'baseline'),'Color','b','LineWidth',2)
line(w,feval(fit,w,'eOPN'),'Color','g','LineWidth',2)
title('Fitted Regression Lines by Conditions')
anova(fit)

%% Regression for eOPN M2 PGD to Neural Trajectory
load('Y:\Hammad\Ephys\LeverTask\eOPNM2\23949eOPNRbp4\Day4\baselineRTSpeed.mat')
DepPGD = dPGD;
load('Y:\Hammad\Ephys\LeverTask\eOPNM2\23949eOPNRbp4\Day4\eOPNRTSpeed.mat')
DepPGD = [DepPGD;dPGD];

experiment = table(DepPGD,IndNeuralSpeed,Condition); % Make a table of how reactiontime depends on neural speed
experiment.Condition = categorical(experiment.Condition);

fit = fitlm(experiment,'DepPGD~IndNeuralSpeed*Condition')

w = linspace(min(IndNeuralSpeed),max(IndNeuralSpeed));
figure()
gscatter(IndNeuralSpeed,DepPGD,Condition,'bgr','x.o')
line(w,feval(fit,w,'baseline'),'Color','b','LineWidth',2)
line(w,feval(fit,w,'eOPN'),'Color','g','LineWidth',2)
title('Fitted Regression Lines by Conditions')
anova(fit)

%% Regression for eOPN in Thalamus

load('Y:\Hammad\Ephys\LeverTask\eOPNM2\23949eOPNRbp4\Day6\baselineRTSpeed.mat')
Condition = repmat("baseline",length(PQ),1); % Make the categorical var
IndNeuralSpeed = PQ; % Investigate the effect of NT Speed on RT
DepReactionT  = rawreactionTime; % Output (RT)
load('Y:\Hammad\Ephys\LeverTask\eOPNM2\23949eOPNRbp4\Day6\eOPNRTSpeed.mat')
Condition = [Condition;repmat("eOPN",length(PQ),1)]; % Make the categorical var
IndNeuralSpeed = [IndNeuralSpeed;PQ]; % Investigate the effect of NT Speed on RT
DepReactionT  = [DepReactionT;rawreactionTime]; % Output (RT)
experiment = table(DepReactionT,IndNeuralSpeed,Condition); % Make a table of how reactiontime depends on neural speed
experiment.Condition = categorical(experiment.Condition);

fit = fitlm(experiment,'DepReactionT~IndNeuralSpeed*Condition')

w = linspace(min(IndNeuralSpeed),max(IndNeuralSpeed));
figure()
gscatter(IndNeuralSpeed,DepReactionT,Condition,'bgr','x.o')
line(w,feval(fit,w,'baseline'),'Color','b','LineWidth',2)
line(w,feval(fit,w,'eOPN'),'Color','g','LineWidth',2)
title('Fitted Regression Lines by Conditions')
anova(fit)

figure()
scatter(IndNeuralSpeed(1:87),DepReactionT(1:87),'k','filled'),set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
ylabel('Reaction time (s)'),xlabel('Neural trajectory speed')
figure()
scatter(IndNeuralSpeed(88:end),DepReactionT(88:end),'o','filled'),set(gca,'TickDir','out'),set(gca,'fontsize',16),box off
ylabel('Reaction time (s)'),xlabel('Neural trajectory speed')

fit = fitlm(IndNeuralSpeed(1:87),DepReactionT(1:87))
fit = fitlm(IndNeuralSpeed(88:end),DepReactionT(88:end))
%% Example script
load carsmall
figure()
gscatter(Weight,MPG,Model_Year,'bgr','x.o')
title('MPG vs. Weight, Grouped by Model Year')
cars = table(MPG,Weight,Model_Year);
cars.Model_Year = categorical(cars.Model_Year);
fit = fitlm(cars,'MPG~Weight*Model')
anova(fit)
w = linspace(min(Weight),max(Weight));
figure()
gscatter(Weight,MPG,Model_Year,'bgr','x.o')
line(w,feval(fit,w,'70'),'Color','b','LineWidth',2)
line(w,feval(fit,w,'76'),'Color','g','LineWidth',2)
line(w,feval(fit,w,'82'),'Color','r','LineWidth',2)
title('Fitted Regression Lines by Model Year')




