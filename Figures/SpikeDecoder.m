%% Load CNN spike decoder data in
filename = 'C:\Users\khan332\Box\ET_NNT\Results\spk_total_confidences.npy';
data = readNPY(filename);
%% Plot out data
t = -1.5:0.001:1.5;
figure;
dat = squeeze(data(:,:,1));
dat(dat(:,2000)<0.01,:) = [];
subplot(131),plot(t,dat','k'); %Hits
hold on, plot(t,mean(dat,1),'r','linewidth',3)
axis square, box off, set(gca,'tickdir','out','fontsize',12),ylabel('Confidence'),xlim([-0.5 1.5])

dat = squeeze(data(:,:,2));
dat(dat(:,2000)<0.01,:) = [];
subplot(132),plot(t,dat','k'); %FA
hold on, plot(t,mean(dat,1),'r','linewidth',3)
axis square, box off, set(gca,'tickdir','out','fontsize',12),xlabel('Time'),xlim([-0.5 1.5])

dat = squeeze(data(:,:,3));
dat(dat(:,2000)<0.01,:) = [];

subplot(133),plot(t,dat','k'); %Miss
hold on, plot(t,mean(dat,1),'r','linewidth',3)
axis square, box off, set(gca,'tickdir','out','fontsize',12),xlim([-0.5 1.5])
yline(0.3,'k')