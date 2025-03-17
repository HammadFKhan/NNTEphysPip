%% Example lever cue session from animal 
clear
load('Y:\Hammad\Ephys\LeverTask\Cue\notagSOM_M2Grid_M1Shank\Day4\IntanBehaviourRT.mat')
%%
dt = mean(diff(IntanBehaviour.time));
time = dt:dt:length(IntanBehaviour.time)*dt;
f = figure(1)
clf
subplot(311),plot(IntanBehaviour.rewardTrace),hold on
axis off
subplot(312),plot(smoothdata(IntanBehaviour.leverTrace,'gaussian',1000)),hold on
yline(0.1,'r'),box off
subplot(313),plot(IntanBehaviour.cueTrace)
axis off
linkaxes
xlim([1.133*10^6 1.21*10^6 ])
ylim([0 0.25])
box off
f.Position = [681 559 860 420];