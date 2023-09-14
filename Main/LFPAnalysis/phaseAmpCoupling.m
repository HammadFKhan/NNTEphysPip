%% Phase amplitude coupling 
% Analysis for time invariate phase amplutude coupling
% data should be segmented for one electrode for specific behaviour trial
% Code reference from mikexcohen.com
%%------ INPUT ------------ %%
% data: single array of data 1xsamples
% srate : sampling rate in Hz
% phase_freq: array of phase frequencies to look at
% ampl_freq: array of amplitude frequencies to look at
function phaseamp = phaseAmpCoupling(data,srate,phas_freqs,ampl_freqs)
showplot = 0; %plot output

if ~exist('phas_freqs','var')
    % define frequencies for phase and for amplitude
    phas_freqs =  5:2:30;
    disp('Using preset phase freq.')
end
if ~exist('ampl_freqs','var')
    ampl_freqs = 40:5:100;
    disp('Using preset amplitude freq.')
end

% number of iterations used for permutation testing
n_iter = 200; 
npnts = size(data,2);
% initialize output phase-amplitude matrix
phaseamp = zeros(length(phas_freqs),length(ampl_freqs));


% loop over frequencies for phase
fprintf('Calculating PAC...')
for lower_fi=1:length(phas_freqs)
    
    % get phase values
    phasefilt = filterFGx(data,srate,phas_freqs(lower_fi),phas_freqs(lower_fi)*.4);
    phase = angle(hilbert(phasefilt));
    
    for upper_fi=1:length(ampl_freqs)
        
        % get power values (note: 'power' is a built-in function so we'll name this variable 'amp')
        ampfilt = filterFGx(data,srate,ampl_freqs(upper_fi),ampl_freqs(upper_fi)*.78);
        amplit = abs(hilbert(ampfilt)).^2;
        
        % calculate observed modulation index
        modidx = abs(mean(amplit.*exp(1i*phase)));

        % now use permutation testing to get Z-value
        bm = zeros(1,length(n_iter));
        for bi=1:n_iter
            cutpoint = randsample(round(npnts/10):round(npnts*.9),1);
            bm(bi) = abs(mean(amplit([ cutpoint:end 1:cutpoint-1 ]).*exp(1i*phase)));
        end

        % the value we use is the normalized distance away from the mean of
        % boot-strapped values
        phaseamp(lower_fi,upper_fi) = (modidx-mean(bm))/std(bm);
    end % end upper frequency loop (for amplitude)
end % end lower frequency loop (for phase)
fprintf('done\n')
if showplot
    % plot it! let's try a contour map
    figure(1), clf
    subplot(3,1,[1 2])
    contourf(phas_freqs,ampl_freqs,phaseamp',40,'linecolor','none')
    set(gca,'clim',[-3 3]),colormap(jet)
    xlabel('Frequency for phase')
    ylabel('frequency for amplitude')
    title('Map of phase-amplitude coupling')
    subplot(313),plot(data,'k','LineWidth',2),box off,set(gca,'tickdir','out')
end
% note: red colors mean more coupling