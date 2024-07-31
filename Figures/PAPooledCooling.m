%% LFP PA Phase for Cooling
% Pooled responses across animals
files = dir(fullfile('D:\M1Cooling\LFP\','*.mat'));
for fileNum = 1:length(files)
    fName = fullfile(files(fileNum).folder,files(fileNum).name);
    disp(['Loading ' fName '...'])
    load(fName)
    % Extract PA structure/values
    PAall(fileNum).CooledPA = LFP.probe1.CooledPA;
    PAall(fileNum).BaselinePA = LFP.probe1.BaselinePA;
end
%% Main analysis here after completing loop

% Calculate statistical PA that is depth specific and as a function of
% trials

PAhitL23 = [];PAhitL5 = []; PAmissL23 = [];  PAmissL5 = [];PAMIFAL23 = [];PAMIFAL5 = [];
for n = 1:length(PAall)
    PAall(n).PA = PAall(n).BaselinePA;
    dat = PAall(n).PA.Hit;    
    sz = size(dat);
    lr = ceil(sz(1)/2);
    
    PAhitL23(n,:) = mean(squeeze(dat(1:lr-1,:,:)));
    PAhitL5(n,:) = mean(squeeze(dat(lr:end,:,:)));
    
    dat = PAall(n).PA.Miss;
    PAmissL23(n,:) = mean(squeeze(dat(1:lr-1,:,:)));
    PAmissL5(n,:) = mean(squeeze(dat(lr:end,:,:)));
    
    dat = PAall(n).PA.MIFA;
    PAMIFAL23(n,:) = mean(squeeze(dat(1:lr-1,:,:)));
    PAMIFAL5(n,:) = mean(squeeze(dat(lr:end,:,:)));
end
%%
time = (-1500:1500)/1000;

dat = smoothdata([PAhitL23;PAhitL5],2,'movmean',75);
figure,plot(time,mean(dat),'k','Linewidth',2),hold on
plot(time,mean(dat)+std(dat)/sqrt(size(dat,1)),'k')
plot(time,mean(dat)-std(dat)/sqrt(size(dat,1)),'k')
box off,set(gca,'tickdir','out','fontsize',14)
xlabel('Time from cue(s)'),ylabel('Phase alignment'),axis square
xlim([-1.5 1.5]),ylim([0.0 0.35])

dat = smoothdata([PAmissL23;PAmissL5],2,'movmean',75);
figure,plot(time,mean(dat),'k','Linewidth',2),hold on
plot(time,mean(dat)+std(dat)/sqrt(size(dat,1)),'k')
plot(time,mean(dat)-std(dat)/sqrt(size(dat,1)),'k')
box off,set(gca,'tickdir','out','fontsize',14)
xlabel('Time from cue(s)'),ylabel('Phase alignment'),axis square
xlim([-1.5 1.5]),ylim([0.0 0.35])       

dat = smoothdata([PAMIFAL23;PAMIFAL5],2,'movmean',75)/2;
figure,plot(time,mean(dat),'k','Linewidth',2),hold on
plot(time,mean(dat)+std(dat)/sqrt(size(dat,1)),'k')
plot(time,mean(dat)-std(dat)/sqrt(size(dat,1)),'k')
box off,set(gca,'tickdir','out','fontsize',14)
xlabel('Time from cue(s)'),ylabel('Phase alignment'),axis square
xlim([-1.5 1.5]),ylim([0.0 0.35])
%%  Evoked pre and post response
%% Layer Specific box plots
%max cued evoked response across layers
PAhitMax = [max(PAhitL23(:,1500:end),[],2),max(PAhitL5(:,1500:end),[],2);];
figure,customBoxplot(PAhitMax),ylim([0 1])
box off,set(gca,'tickdir','out','fontsize',14)
ylabel('Phase alignment of Layers'),axis square

PAmissMax = [max(PAmissL23(:,1500:end),[],2),max(PAmissL5(:,1500:end),[],2);];
figure,customBoxplot(PAmissMax),ylim([0 1])
box off,set(gca,'tickdir','out','fontsize',14)
ylabel('Phase alignment of Layers'),axis square

PAMIFAMax = [max(PAMIFAL23(:,1500:end),[],2),max(PAMIFAL5(:,1500:end),[],2);];
figure,customBoxplot(PAMIFAMax),ylim([0 1])
box off,set(gca,'tickdir','out','fontsize',14)
ylabel('Phase alignment of Layers'),axis square