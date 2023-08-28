out=import_wcp;
Vm = out.S{1};
Im = out.S{2};
T = out.T;
dt = out.T(2)-out.T(1);% time
Fs = 1/dt; 
%% Vm
% Fix DC offset
VmTotal = [];
% Vm = Vm-min(Vm,[],'all');
% Vm = Vm-65;
for i = 800:size(Vm,2)
    VmTotal = [VmTotal;Vm(:,i)];
end

% VmTotal = VmTotal-min(VmTotal,[],'all');
% VmTotal = VmTotal-62;
figure,plot(0:dt:(length(VmTotal)-1)/Fs,VmTotal);xlim([0 length(VmTotal)/Fs])
% figure,plot(0:dt:(length(VmfiltBeta)-1)/Fs,VmfiltBeta);xlim([0 length(VmTotal)/Fs])
%% Im
% Fix DC offset
ImTotal = [];
Im = Im-min(Im,[],'all');
for i = 1:size(Im,2)
    ImTotal = [ImTotal;Im(:,i)];
end

figure,plot(0:dt:(length(ImTotal)-1)/Fs,ImTotal);xlim([0 length(ImTotal)/Fs])

%% STA
thresh = find(VmTotal>83);
tri = diff(thresh);
trig = find(tri~=1);
COM = thresh(trig+1,1);
for i = 1:length(trig)
    win(:,i) = VmTotal(COM(i)-(0.001*Fs):COM(i)+(0.002*Fs),1);
end
% upWin = interp1(1:size(win,1),win,1:0.05:size(win,1),'spline');
%% Delete action potentials for subthreshold
subThreshold = VmTotal;
%%

thresh = find(subThreshold>-130);
tri = diff(thresh);
trig = find(tri~=1);
COM = thresh(trig,1);

subThreshold(COM(2:9000)-(0.001*Fs):COM(2:9000)+(0.002*Fs)) = [];



figure,plot(0:dt:(length(subThreshold)-1)/Fs,subThreshold,'r');xlim([0 length(subThreshold)/Fs]),box off

%%
figure,plot(upWin(:,1:150)),hold on
plot(mean(upWin(:,1:150),2),'k','LineWidth',3),axis off
%% Subthreshold Dynamics
% subthresh = Intracellular2<-140;
subthresh = Intracellular2(subthresh);
Fc = [1 100];
Wn = Fc./(Fs/2);
b = fir1(10000,Wn,'bandpass');
subthreshold = filtfilt(b,1,LFPtest);
%% Filter
filtered_dataIntra = customFilt(subThreshold',Fs,[12 30]);
%%
Intra = IntrabetaBurstDetection(filtered_dataIntra',Fs);
%%
figure
for i = 1:81
    subplot(9,9,i),plot(Intra.betaTrace{i}),axis off
end
%% Wavelet 
[wavelet, f] = cwt(Vm(:,1),Fs,'FrequencyLimit',[10 100]);
figure,imagesc(0:dt:(length(subThreshold)-1)/Fs,f,abs(wavelet));colormap(jet);axis xy
%%
[peakAlign,norm,f,stats] = IntrabetaAnalysis(Intra);
%%
figure,plot(filtered_dataIntra),hold on
for i = 1:87
    xline(LFP.betaBurst.detectedBeta(i,2)*Fs,'r--','LineWidth',1)
end
%%
figure,
subplot(2,1,1),plot(0:dt:(length(VmTotal)-1)/Fs,VmTotal);xlim([0 length(VmTotal)/Fs])
subplot(2,1,2),plot(0:dt:(length(VmTotal)-1)/Fs,filtered_dataIntra);xlim([0 length(VmTotal)/Fs]),linkaxes
% subplot(2,1,2),plot(subthresh);xlim([0 length(VmTotal)/Fs])
%%
