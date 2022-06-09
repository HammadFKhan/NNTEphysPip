out=import_wcp;
Vm = out.S{1};
Im = out.S{2};
T = out.T;
dt = out.T(2)-out.T(1);% time
Fs = 1/dt; 
%% Vm
% Fix DC offset
VmTotal = [];
Vm = Vm-min(Vm,[],'all');
Vm = Vm-65;
for i = 1:size(Vm,2)
    VmTotal = [VmTotal;Vm(:,i)];
end

figure,plot(0:dt:(length(VmTotal)-1)/Fs,VmTotal);xlim([0 length(VmTotal)/Fs])
%% Im
% Fix DC offset
ImTotal = [];
Im = Im-min(Im,[],'all');
for i = 1:size(Im,2)
    ImTotal = [ImTotal;Im(:,i)];
end

figure,plot(0:dt:(length(ImTotal)-1)/Fs,ImTotal);xlim([0 length(ImTotal)/Fs])
%% STA
thresh = find(VmTotal>10);
tri = diff(thresh);
trig = find(tri~=1);
COM = thresh(trig+1,1);
for i = 1:length(trig)
    win(:,i) = VmTotal(COM(i)-6:COM(i)+15,1);
end
upWin = interp1(1:size(win,1),win,1:0.05:size(win,1),'spline');
%% Delete action potentials for subthreshold
subThreshold = VmTotal;
for i = 1:20
    thresh = find(subThreshold>0);
    tri = diff(thresh);
    trig = find(tri~=1);
    COM = thresh(trig,1);
    subThreshold(COM(1)-6:COM(1)+15) = [];
end
%%
figure,plot(upWin(:,1:20)),hold on
plot(mean(upWin(:,1:20),2),'k','LineWidth',3),axis off
%% Subthreshold Dynamics
subthresh = VmTotal<30;
subthresh = VmTotal(subthresh);
Fc = [1 15];
Wn = Fc./(Fs/2);
b = fir1(50,Wn,'bandpass');
subThreshold = filtfilt(b,1,subthresh);
%% Filter
filtered_data = customFilt(subThreshold',Fs,[10 30]);
%%
LFP = IntrabetaBurstDetection(filtered_data',Fs);
%%
figure
for i = 1
    subplot(1,1,i),plot(LFP.betaTrace{i}),axis off
end
%% Wavelet 
[wavelet, f] = cwt(filtered_data,Fs,'FrequencyLimit',[10 30]);

%%
figure,
subplot(2,1,1),plot(0:dt:(length(VmTotal)-1)/Fs,VmTotal);xlim([0 length(VmTotal)/Fs])
subplot(2,1,2),plot(0:dt:(length(VmTotal)-1)/Fs,filtered_data);xlim([0 length(VmTotal)/Fs]),linkaxes
% subplot(2,1,2),plot(subthresh);xlim([0 length(VmTotal)/Fs])