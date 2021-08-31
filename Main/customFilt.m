function customFilt(data,Fs,~)
%% band-pass filtering
% specify Nyquist freuqency
nyquist = Fs/2;
nbchan = size(data,1);
pnts = size(data,2);
times = 0:(pnts/Fs);
% filter frequency band
filtbound = [4 10]; % Hz

% transition width
trans_width = 0.2; % fraction of 1, thus 20%

% filter order
filt_order = round(3*(Fs/filtbound(1)));

% frequency vector (as fraction of Nyquist
ffrequencies  = [ 0 (1-trans_width)*filtbound(1) filtbound (1+trans_width)*filtbound(2) nyquist ]/nyquist;

% shape of filter (must be the same number of elements as frequency vector
idealresponse = [ 0 0 1 1 0 0 ];

% get filter weights
filterweights = firls(filt_order,ffrequencies,idealresponse);

% plot for visual inspection
figure(), clf
subplot(211)
plot(ffrequencies*nyquist,idealresponse,'k--o','markerface','m')
set(gca,'ylim',[-.1 1.1],'xlim',[-2 nyquist+2])
xlabel('Frequencies (Hz)'), ylabel('Response amplitude')

subplot(212)
plot((0:filt_order)*(1000/Fs),filterweights)
xlabel('Time (ms)'), ylabel('Amplitude')

% apply filter to data
filtered_data = zeros(nbchan,pnts);
for chani=1:nbchan
    filtered_data(chani,:) = filtfilt(filterweights,1,double(data(chani,:)));
end

figure(), clf
plot(times,squeeze(data(3,:)))
hold on
plot(times,squeeze(filtered_data(3,:)),'r','linew',2)
xlabel('Time (ms)'), ylabel('Voltage (\muV)')
legend({'raw data';'filtered'})