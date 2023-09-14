function info = tapered_spike_coherence_modified(spikea,spikeb,interval,window,showplot)
%********* JUDE MITCHELL (jude@salk.edu):  7/18/2009
%********* HAMMAD KHAN (khan332@purdue.edu): 7/23/2023
%******** function info = tapered_spike_coherence(data,data2,interval,window,showplot)
%******* example: info =
%                   spikeCoherence(trial).spikecoherence = ...
%            tapered_spike_coherence_modified(spikeCoherence(trial).spikea,spikeCoherence(trial).spikeb,[],800,0);
%
%**********
%***** general: computes the coherence using Multiple tapers
%** inputs:
%********** data - all relevent data fields for a neuron
%******     spikea - cell structure of trials with discretized AP 
%****** 
%****** 
%****** 
%******                
%******
%********** spikeb - data for second unit
%********** interval - interval of analysis, 1xN array of timepoints
%********** window - duration of LFP window around each spike for FFT (start with 800)
%********** showplot - to plot out results
%*****
%****** showplot - set this to 1 to see the results plotted out
%*****
%******  output:
%******      info.acoho  - attended coherence [1 x F]
%******      info.ascoho - std from permutes of estimate
%******      info.arcoho  - shuffled trials, attended coherence [1 x F]
%******      info.arscoho - shuffled trials, std from permutes of estimate
%******      info.afreq - frequencies [1 x F]
%******      info.ucoho  - attended coherence [1 x F]
%******      info.uscoho - std from jacknife of estimate
%******      info.urcoho  - shuffled trials, attended coherence [1 x F]
%******      info.usrcoho - shuffled trials, std from jacknife
%******      info.ufreq - frequencies [1 x F]

%datapath = 'export_datapair';


%************ plot the spike rasters and mean firing rates to sanity check
%************ and also show the local fields and mean local fields
if (showplot == 1)
    % Make some nice spike plots of data
    spiker = spikea;
    
    figure;
    %****************
    disp('Plotting rastergrams (slow) ...');
    subplot('position',[0.1 0.4 0.4 0.55]);
    make_nice_spike_raster(spiker);
    V = axis;
    axis([0 size(spikea{1},2) V(3) V(4)]);
    grid on;
    ylabel('Trial Number');
    title(fprintf('Unit rasters'));
    %*************
    subplot('position',[0.1 0.1 0.4 0.3]);
    smooth_window = 25;  % give sigma of 12.5ms
    make_nice_mean_raster(spiker,smooth_window);
    V = axis;
    axis([0 size(spikea{1},2) V(3) V(4)]);
    plot([interval(1),interval(1)],[V(3),V(4)],'k-'); hold on;
    plot([interval(end),interval(end)],[V(3),V(4)],'k-'); hold on;
    ylabel('Firing Rate');
    xlabel('Time (ms)');
    
    spiker = spikeb;
    %****************
    disp('Plotting rasters2 per trial (slow) ...');
    subplot('position',[0.57 0.4 0.4 0.55]);
    make_nice_spike_raster(spiker);
    V = axis;
    axis([0 size(spikeb{1},2) V(3) V(4)]);
    grid on;
    ylabel('Trial Number');
    title(fprintf('Unit rasters'));
    %*************
    subplot('position',[0.57 0.1 0.4 0.3]);
    smooth_window = 25;  % give sigma of 12.5ms
    make_nice_mean_raster(spiker,smooth_window);
    V = axis;
    axis([0 size(spikeb{1},2) V(3) V(4)]);
    plot([interval(1),interval(1)],[V(3),V(4)],'k-'); hold on;
    plot([interval(end),interval(end)],[V(3),V(4)],'k-'); hold on;
    ylabel('Firing Rate');
    xlabel('Time (ms)');
    
end


%**************** now call for the coherence estimates ***********
Win = window;  % use specified window size
%************************************************
M = 20;  % repeat analysis several times doing a
% and random repeat of the rate matching
% between the attention conditions
%******** find min number of spikes per condition ******
minspikecnta = 10000000000000;
minspikecntb = 10000000000000;
CNUM = length(spikea); % Number of experimental conditions or trials
for cubo = 1:CNUM
    %***************************
    it = sum( sum( spikea{cubo}(:,interval) ) );
    if (it < minspikecnta)
        minspikecnta = it;
    end
    %***************************
    it = sum( sum( spikeb{cubo}(:,interval) ) );
    if (it < minspikecntb)
        minspikecntb = it;
    end
    %****************************
end
%****** finds the minimum spike count in any
%****** attention condition for each unit, then
%****** it will match counts across attention conditions
%****************************
zcoho = cell(1,CNUM);
zscoho = cell(1,CNUM);
zphaso = cell(1,CNUM);
zsphaso = cell(1,CNUM);
zspika_pow = cell(1,CNUM);
zsspika_pow = cell(1,CNUM);
zspike_pow = cell(1,CNUM);
zsspike_pow = cell(1,CNUM);
zrcoho = cell(1,CNUM);
zrphaso = cell(1,CNUM);
zfreq = cell(1,CNUM);
%**********************************************
for cubo = 1:CNUM
    fprintf('Computing coherence on group %d trials\n',cubo);
    for ii = 1:M
        fprintf('Computing interation %d...', ii)
        spikeholea = spikea{cubo}(:,interval);
        spikeholeb = spikeb{cubo}(:,interval);
        
        if (1)  % do rate matching
            %**** match the spike rate to minspikecnt ******
            y = find( spikeholea > 0 );
            yn = size(y,1);  % number of spikes
            permy = randperm(yn);
            tossout = yn - minspikecnta;
            spikeholea( y(permy(1:tossout)) ) = 0;
            %**** match the spike rate for other spike train ***
            y = find( spikeholeb > 0 );
            yn = size(y,1);  % number of spikes
            permy = randperm(yn);
            tossout = yn - minspikecntb;
            spikeholeb( y(permy(1:tossout)) ) = 0;
        end
        %************ compute coherence from resampled version
        [icoho,iscoho,iphaso,isphaso,ispika_pow,isspika_pow,...
            ispike_pow,isspike_pow,ircoho,irphaso,ifreq] = ...
            taper_coherence(spikeholea,spikeholeb,Win);
        %*********** pool the results *******************
        zcoho{cubo} = [ zcoho{cubo} ; icoho];
        zscoho{cubo} = [ zscoho{cubo} ; iscoho];
        zphaso{cubo} = [ zphaso{cubo} ; iphaso];
        zsphaso{cubo} = [ zsphaso{cubo} ; isphaso];
        zspika_pow{cubo} = [ zspika_pow{cubo} ; ispika_pow];
        zsspika_pow{cubo} = [ zsspika_pow{cubo} ; isspika_pow];
        zspike_pow{cubo} = [ zspike_pow{cubo} ; ispike_pow];
        zsspike_pow{cubo} = [ zsspike_pow{cubo} ; isspike_pow];
        zrcoho{cubo} = [ zrcoho{cubo} ; ircoho];
        zrphaso{cubo} = [ zrphaso{cubo} ; irphaso];
        zfreq{cubo} = [ zfreq{cubo} ; ifreq];
        fprintf('done!\n')
        
    end
    if ii == 1
        coho{cubo} = zcoho{cubo};
        scoho{cubo} = zscoho{cubo};
        phaso{cubo} = zphaso{cubo};
        sphaso{cubo} = zsphaso{cubo};
        spika_pow{cubo} = zspika_pow{cubo};
        sspika_pow{cubo} = zsspika_pow{cubo};
        spike_pow{cubo} = zspike_pow{cubo};
        sspike_pow{cubo} = zsspike_pow{cubo};
        rcoho{cubo} = zrcoho{cubo};
        rphaso{cubo} = zrphaso{cubo};
        freq{cubo} = zfreq{cubo};
    else
        coho{cubo} = mean( zcoho{cubo} );
        scoho{cubo} = mean( zscoho{cubo} );
        phaso{cubo} = mean( zphaso{cubo} );
        sphaso{cubo} = mean( zsphaso{cubo} );
        spika_pow{cubo} = mean( zspika_pow{cubo} );
        sspika_pow{cubo} = mean( zsspika_pow{cubo} );
        spike_pow{cubo} = mean( zspike_pow{cubo} );
        sspike_pow{cubo} = mean( zsspike_pow{cubo} );
        rcoho{cubo} = mean( zrcoho{cubo} );
        rphaso{cubo} = mean( zrphaso{cubo} );
        freq{cubo} = mean( zfreq{cubo} );
    end
    
end
%*********** store results *********
info.coho = coho;
info.scoho = scoho;
info.phaso = phaso;
info.sphaso = sphaso;
info.spika_pow = spika_pow;
info.sspika_pow = sspika_pow;
info.spike_pow = spike_pow;
info.sspike_pow = sspike_pow;
info.rcoho = rcoho;
info.rphaso = rphaso;
info.freq = freq;
%***********************************

if (showplot == 1)
    figure;
    colo = 'rbgy';
    
    subplot(2,2,1);
    for ii = 1:CNUM
        H = semilogx(freq{ii},coho{ii},[colo(ii),'-']); hold on;
        set(H,'Linewidth',2);
        H = semilogx(freq{ii},(coho{ii}+scoho{ii}),[colo(ii),':']); hold on;
        H = semilogx(freq{ii},(coho{ii}-scoho{ii}),[colo(ii),':']);
        H = semilogx(freq{ii},rcoho{ii},[colo(ii),'--']); hold on;
        set(H,'Linewidth',1);
    end
    ylabel('Coherence Magnitude');
    xlabel('Frequency (Hz)');
    title(sprintf('Spike A with Spike B'));
    
    subplot(2,2,2);
    for ii = 1:CNUM
        H = plot(freq{ii},phaso{ii},[colo(ii),'-']); hold on;
        set(H,'Linewidth',2);
        H = plot(freq{ii},(phaso{ii}+sphaso{ii}),[colo(ii),':']); hold on;
        H = plot(freq{ii},(phaso{ii}-sphaso{ii}),[colo(ii),':']);
        %H = plot(freq{ii},rphaso{ii},[colo(ii),'--']); hold on;
        set(H,'Linewidth',1);
    end
    ylabel('Coherence Angle');
    xlabel('Freq');
    title(sprintf('Spike A with Spike B'));
    
    subplot(2,2,3);
    for ii = 1:CNUM
        H = plot(freq{ii},spika_pow{ii},[colo(ii),'-']); hold on;
        set(H,'Linewidth',2);
        H = plot(freq{ii},(spika_pow{ii}+sspika_pow{ii}),[colo(ii),':']); hold on;
        H = plot(freq{ii},(spika_pow{ii}-sspika_pow{ii}),[colo(ii),':']);
    end
    ylabel('Spike Unit A Pow');
    xlabel('Frequency (Hz)');
    title(sprintf('Spike A with Spike B'));
    
    subplot(2,2,4);
    for ii = 1:CNUM
        H = plot(freq{ii},spike_pow{ii},[colo(ii),'-']); hold on;
        set(H,'Linewidth',2);
        H = plot(freq{ii},(spike_pow{ii}+sspike_pow{ii}),[colo(ii),':']); hold on;
        H = plot(freq{ii},(spike_pow{ii}-sspike_pow{ii}),[colo(ii),':']);
    end
    ylabel('Spike Unit B Pow');
    xlabel('Frequency (Hz)');
    title(sprintf('Spike B with Spike A'));
    
end

return;


function [coho,scoho,phaso,sphaso,spika_pow,sspika_pow,spike_pow,sspike_pow,...
    rcoho,rphaso,ifcoher] = taper_coherence(spika,spiko,Win)
%******* [coho,phaso,spika_pow,spike_pow,...
%         rcoho,rphaso,ifcoher] = taper_coherence(speco,spiko,Win)
%****** computes the coherence value using taper estimates by breaking
%****** up data intervals into overlapping intervals of size Win
%****** inputs:
%******    spika - a MxT array where M is trials and T is time in ms
%******    spiko - a MxT array with same dimensions, 0 or 1 for spikes
%******    Win - a window around which to compute LFP segments
%****** outputs:
%******    coho - 1xF array of magnitude coherence per frequencies F
%******    phaso - 1xF array of phase values
%******    spika_pow - 1xF spike unit a power spectra
%******    spike_pow - 1xF spike power spectra (normalized via rate)
%******    ifcoher - list of frequency values, 1xF

%*********** prep multi-taper fft ***************
TT = size(spika,2);
N = Win;
Trials = size(spika,1);
%*******************
TTa = 1+floor(Win/2);   %tighten up interval of analysis so Win does not fall outside
TTb = TT-floor(Win/2);
%******** tapers with smoothing
W = 2.5;  % 2.5 hz
NW = floor( 2 * ((Win/1000)*W))/2;
tap = dpss(Win,NW);
tapers = tap';
KW = floor( 2*NW - 1);
fprintf('Computing with %d tapers\n',KW);
%********** compute range of frequencies for FFT analysis
pad = 0;
fpass = [0.001 0.088];   % our LFP has a restricted range due to hardware filter
Fs = 1;
nfft=2^(nextpow2(N));
df=Fs/nfft;
freqreal=0:df:(Fs-df); % all possible frequencies
findx=find(freqreal >= fpass(1) & freqreal <= fpass(end));
f=freqreal(findx);
%**********************

for repo = 1:2
    zcross_pow = [];
    zispika_pow = [];
    zispike_pow = [];
    zspike_count = [];
    zspika_count = [];
    %********************************
    %******* run through all segments, allow overlap via Win/4
    if (repo == 1)  % repo 1, use aligned trials
        tr_perm = 1:Trials;
    else            % repo 2, do same procedure on shifted trials
        mid = floor(Trials/8) + floor( rand * (3*Trials/4) );
        tr_perm = [mid:Trials,1:(mid-1)];  % shift-predictor
        % tr_perm = randperm(Trials);
    end
    for tr = 1:Trials
        %***** get sum from each trial, use to make jacknifes
        cross_pow = [];
        ispika_pow = [];
        ispike_pow = [];
        spike_count = [];
        spika_count = [];
        %***********************
        TTa = 1;
        TTb = Win;
        Tstep = floor(Win/4);
        while (TTb <= TT)
            %*************
            aspiker = spika(tr,TTa:TTb);
            naspiker = aspiker - mean(aspiker);  % set DC to zero
            nspiker = spiko(tr_perm(tr),TTa:TTb);
            spiker = nspiker - mean(nspiker); % set DC to zero
            %**************
            for kt = 1:KW  % repeat over different tapers
                aspikero = naspiker .* tapers(kt,:);
                spikero = spiker .* tapers(kt,:);
                
                J1 = fft(aspikero,nfft);
                J2 = fft(spikero,nfft);
                
                if (isempty(ispika_pow))  % first data point set it
                    ispika_pow = J1 .* conj(J1);  % power unit A
                    ispike_pow = J2 .* conj(J2);  % power unit B
                    cross_pow = J1 .* conj(J2);   % cross-correlation A and B
                    spike_count = sum(nspiker);
                    aspike_count = sum(aspiker);
                else                      % remaining points add it up
                    ispika_pow = ispika_pow + (J1 .* conj(J1));
                    ispike_pow = ispike_pow + (J2 .* conj(J2));
                    cross_pow = cross_pow + (J1 .* conj(J2));
                    spike_count = spike_count + sum(nspiker);
                    aspike_count = aspike_count + sum(aspiker);
                end
            end
            %***************
            TTa = TTa + Tstep;
            TTb = TTb + Tstep;
        end
        %************** store estimate from each trial ***********
        %************** so we can use jack-knife error bars ******
        zcross_pow = [ zcross_pow ; cross_pow ];
        zispika_pow = [ zispika_pow ; ispika_pow ];
        zispike_pow = [ zispike_pow ; ispike_pow];
        zspike_count = [ zspike_count ; spike_count];
        zspika_count = [ zspika_count ; aspike_count];
    end
    
    if (repo == 1)
        
        %******* jacknife the coherence estimate for error bars ***
        cohjack = [];
        phajack = [];
        spajack = [];
        spijack = [];
        %**************
        for tr = 1:Trials
            y = [1:(tr-1),(tr+1):Trials];
            cross_pow = sum( zcross_pow(y,:) );
            ispike_pow = sum( zispike_pow(y,:) );
            ispika_pow = sum( zispika_pow(y,:) );
            spike_count = sum( zspike_count(y) );
            aspike_count = sum( zspika_count(y) );
            %*********************
            coh = cross_pow ./ (sqrt( ispike_pow ) .* sqrt( ispika_pow ));
            coho = abs(coh(findx));
            phaso = angle(coh(findx));
            spika_pow = ispika_pow(findx) / aspike_count;
            spike_pow = ispike_pow(findx) / spike_count;
            ifcoher = 1000 * freqreal(findx);
            %**********************
            cohjack = [cohjack ; coho];
            phajack = [phajack ; coh];
            spajack = [spajack ; spika_pow];
            spijack = [spijack ; spike_pow];
        end
        %**************************
        coho = mean( cohjack );
        scoho = std( cohjack ) * sqrt( Trials-1 );
        pha = mean( phajack );
        phaso = angle( pha( findx) );   %no jacknife for phase
        phabo = [];
        for iz = 1:size(phajack)
            y = angle( phajack(iz,:) );
            it = find( (y(findx) - phaso) > pi);
            y(findx(it)) = y(findx(it)) - (2*pi);
            it = find( (y(findx) - phaso) < -pi );
            y(findx(it)) = y(findx(it)) + (2*pi);
            phabo = [phaso ; y(findx)];
        end
        sphaso = std( phabo ) * sqrt( Trials-1);
        %**********
        spika_pow = mean( spajack );
        sspika_pow = std( spajack ) * sqrt( Trials-1 );
        spike_pow = mean( spijack );
        sspike_pow = std( spijack ) * sqrt( Trials-1 );
        
    else
        
        cross_pow = sum( zcross_pow );
        ispike_pow = sum( zispike_pow );
        ispika_pow = sum( zispika_pow );
        
        coh = cross_pow ./ (sqrt( ispike_pow ) .* sqrt( ispika_pow ));
        
        rcoho = abs(coh(findx));
        rphaso = angle(coh(findx));
        
    end
    
end

return;



%************** below here are some boring graphics routines to plot things

%****************** function to make a nice raster plot ****************
function make_nice_mean_raster(spmat,smooth_window)
%*********** spmat1 and spmat2 are spike matrices of two conditions you wish to compare
%*********** smooth_window ... gaussian smoothing in millisecs
numconds = size(spmat,2);
if (numconds==2)
    colo = [[1,0,0];[0,0,1]];
else
    colo = [[1,0,0];[0,0,1];[0,1,0];[1,1,0]];
end

for k = 1:numconds
    spud = spmat{1,k};
    numtrials(k) = size(spud,1);
    smorate = gauss_smooth(sum( spud(1:numtrials(k),:))/....
        numtrials(k),smooth_window)*1000;
    H = plot(smorate,'k-'); hold on;
    set(H,'Color',colo(k,:));
end

return;


%******************* make a rastergram of the actual spikes on each trial
function make_nice_spike_raster(spmat)

colo = [[1,1,1];[0,0,0]];
colormap(colo);

totspike = [];
for cubo = 1:size(spmat,2)
    totspike = [totspike; (cubo*spmat{1,cubo}) ];
    totspike = [totspike; (cubo*spmat{1,cubo}) ];
    
end
totspike = totspike + 1;
imagesc(totspike);
V = axis;
axis([V(1) V(2) 0 size(totspike,1)]);

return;


%**************************************************************
function output = gauss_smooth(input, window)
% Smoothing function:
% output = smooth(input, window)
% "Window" is the total kernel width.
% Input array must be one-dimensional.

input_dims = ndims(input);
input_size = size(input);
if input_dims > 2 | min(input_size) > 1,
    disp('Input array is too large.');
    return
end

if input_size(2) > input_size(1),
    input = input';
    toggle_dims = 1;
else
    toggle_dims = 0;
end

if window/2 ~= round(window/2),
    window = window + 1;
end
halfwin = window/2;

input_length = length(input);
%********* gauss window +/- 1 sigma
x = -halfwin:1:halfwin;
kernel = exp(-x.^2/(window/2)^2);
kernel = kernel/sum(kernel);

padded(halfwin+1:input_length+halfwin) = input;
padded(1:halfwin) = ones(halfwin, 1)*input(1);
padded(length(padded)+1:length(padded)+halfwin) = ones(halfwin, 1)*input(input_length);

output = conv(padded, kernel);
output = output(window:input_length+window-1);

if toggle_dims == 1,
    output = output';
end

return;



