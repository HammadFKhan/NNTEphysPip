function filtered_data = customFilt(data,Fs,filtbound)
if nargin < 3 || strcmp(filtbound,'')
    error('Include cutoff frequencies [w1 w2]')
end
%% band-pass filtering
% specify Nyquist freuqency
nyquist = Fs/2;
nbchan = size(data,1);
pnts = size(data,2);
times = (1:pnts)/Fs;
% filter frequency band
% filtbound = [4 10]; % Hz

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
figure(1), clf
subplot(211)
plot(ffrequencies*nyquist,idealresponse,'k--o','markerface','m')
set(gca,'ylim',[-.1 1.1],'xlim',[-2 nyquist+2])
xlabel('Frequencies (Hz)'), ylabel('Response amplitude')
xlim([0 2*filtbound(2)])
subplot(212)
plot((0:filt_order)*(1000/Fs),filterweights)
xlabel('Time (s)'), ylabel('Amplitude')

% apply filter to data
filtered_data = zeros(nbchan,pnts);
disp('Applying Filter...')
tic
for chani=1:nbchan
    filtered_data(chani,:) = filtfilt(filterweights,1,data(chani,:));
end
toc

figure(2), clf
plot(times,squeeze(data(1,:)))
hold on
plot(times,squeeze(filtered_data(1,:)),'r','linew',2)
xlabel('Time (s)'), ylabel('Voltage (\muV)')
legend({'raw data';'filtered'})
