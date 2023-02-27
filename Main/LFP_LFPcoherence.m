%% LFP-LFP Coherence

referenceLFP = TimeFreq.tfRun.fullLFP(:,[1:43,45:64],:);
referenceBeta = TimeFreq.tfRun.betaLFP(:,[1:43,45:64],:);
referenceGamma = TimeFreq.tfRun.gammaLFP(:,[1:43,45:64],:);
rightSide = find(s.sorted_probe_wiring(:,2)==0);
rightSide = [rightSide;find(s.sorted_probe_wiring(:,2)==36)];
rightLFP = referenceLFP(:,rightSide,:);
rightBeta = referenceBeta(:,rightSide,:);
rightGamma = referenceGamma(:,rightSide,:);

trial = 1;
figure,stack_plot(rightLFP(:,:,1)',1,2,8192)
figure,stack_plot(rightBeta(:,:,1)',1,2,1024)
figure,stack_plot(rightGamma(:,:,1)',1,2,1024)
%% Coherence of LFP across sites
params.Fs = 1024;
params.tapers = [5 9];
movingwin = [0.5 0.05];
params.pad = 1;
params.fpass = [0 100];
[grad,~]=colorGradient([7 49 97]/255,[110 192 235]/255,25);
figure,hold on
for i = 1:21
data1 = squeeze(rightLFP(:,4,:));
data2 = squeeze(rightLFP(:,15,:));
[C,phi,S12,S1,S2,f]= coherencyc(data1,data2,params);
subplot(2,1,1),plot(f,mean(C,2)),hold on
subplot(2,1,2),plot(f,smoothdata(mean(phi,2),'gaussian',50),'Color',grad(i,:)),title(num2str(i)),ylim([-.2 .2]),hold on
end
set(gca,'TickDir','out')

%% Coherence of LFP across sites
params.Fs = 1024;
params.tapers = [5 9];
movingwin = [0.5 0.05];
params.pad = 1;
params.fpass = [0 100];
[grad,~]=colorGradient([7 49 97]/255,[110 192 235]/255,25);
figure,hold on
data1 = squeeze(rightLFP(:,2,:));
data2 = squeeze(rightLFP(:,15,:));
[C,phi,S12,S1,S2,f]= coherencyc(data1,data2,params);
% subplot(2,1,1),plot(f,mean(C,2)),hold on
lineError(f,smoothdata(mean(phi,2)','gaussian',50),'ste'),title('L23 -> L5'),ylim([-.2 .2]),hold on, box off
set(gca,'TickDir','out')


%% Beta Event
params.Fs = 1024;
params.tapers = [4 7];
movingwin = [0.5 0.05];
params.pad = 1;
params.fpass = [0 100];

t = mLFP{55}([1:43,45:64],:,:);
betaEventRight = t(rightSide,:,:);
figure,hold on
% for i = 1:21
data1 = squeeze(t(2,:,:));
data2 = squeeze(t(4,:,:));
[C,phi,S12,S1,S2,f]= coherencyc(data1,data2,params); % 5 and 15 electrode
% subplot(2,1,1),plot(f,mean(C,2)),hold on
lineError(f,smoothdata(mean(phi,2)','gaussian',10),'ste'),hold on,box off,ylim([-0.2 0.2])
% end
set(gca,'TickDir','out')
%% Sine function
x = 0:pi/25:2*pi;
data1 = sin(x+pi/2); % Compared signal
data2 = sin(x); % Reference signal
data3 = sin(x-pi/2);
figure,plot(data1),hold on
plot(data2);plot(data3);

[C,phi,S12,S1,S2,f]= coherencyc(data3,data2,params); % 5 and 15 electrode
figure,plot(f,smoothdata(mean(phi,2),'gaussian',10),'Color',grad(i,:)),hold on
%%
figure,hold on
for i = 1:21
data1 = squeeze(CSD );
data2 = squeeze(rightLFP(:,15,:));
[C,phi,S12,S1,S2,f]= coherencyc(data1,data2,params);
subplot(2,1,1),plot(f,mean(C,2)),hold on
subplot(2,1,2),plot(f,smoothdata(mean(phi,2),'gaussian',50),'Color',grad(i,:)),hold on
end
set(gca,'TickDir','out')