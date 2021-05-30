function CSD_ripple(ripple_signal)
timestamps = (6.7E5:7E5)';
lfp_frag = ripple_signal(6.7E5:7E5,15:25);
lfp_frag = detrend(lfp_frag')';
lfp_frag = detrend(lfp_frag')';
for ch = 1:size(lfp_frag,2)
    lfp_frag(:,ch) = smooth(lfp_frag(:,ch),0.1,'sgolay');
end
for t = 1:size(lfp_frag,1)
    lfp_frag(t,:) = smooth(lfp_frag(t,:),0.1,'lowess');
end
% calculate CSD 
CSD = diff(lfp_frag,2,2);
% generate output structure
% csd.data = CSD;
% csd.timestamps = timestamps();
% csd.samplingRate = samplingRate;
% csd.channels = channels; 
% csd.params.spat_sm = spat_sm;
% csd.params.temp_sm = temp_sm;
% csd.params.detrend = doDetrend;

figure;
subplot(1,2,1);
cmax = max(max(CSD));
cmin = min(mean(CSD));
contourf(timestamps,1:size(CSD,2),CSD',40,'LineColor','none');hold on;
% imagesc(CSD');hold on;
colormap(jet(10));
caxis([cmin cmax]);
colorbar
xlabel('time (s)');ylabel('channel');title(CSD);

