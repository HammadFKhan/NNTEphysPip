function Spikes = leverPSTH(Spikes,Behaviour)
% Organizes and calculates smoothed firing rates across trials
% This data is organized for coherence measurements
if ~isstruct(Behaviour),error('Behaviour must be structure!');end
if ~isstruct(Spikes), error('Spikes must be a structure!');end

%% Analyze for hit trials
if isfield(Behaviour,'nCueHit')
    fprintf('Analyzing hit trials...\n')
    count = 1;
    trials = {};
    for ii = 1:length(Spikes.Clusters)
        if ~isempty(Spikes.Clusters(ii).spikeTime)
            for i = 1:Behaviour.nCueHit
                st = Behaviour.cueHitTrace(i).LFPtime(1);
                stp = Behaviour.cueHitTrace(i).LFPtime(end);
                temp = zeros(1,length(0:0.001:round((stp-st),3)));
                spiketm = Spikes.Clusters(ii).spikeTime(Spikes.Clusters(ii).spikeTime>=st & Spikes.Clusters(ii).spikeTime<=stp);
                spiked = discretize(spiketm,st:0.001:stp); % 0.001 bin for 1 ms
                if ~isnan(spiked)
                    temp(spiked) = 1;
                end
                trials{count}(i,:) = temp;
                
            end
        else
            trials{count} = [];
        end
        count = count+1;
    end
    %%
    output = make_nice_mean_raster(trials,20,0);
    %%
    Spikes.PSTH.hit.spks = trials;
    Spikes.PSTH.hit.spkRates = output;
end
%% Analyze for miss trials
if isfield(Behaviour,'nCueMiss')
    fprintf('Analyzing miss trials...\n')
    count = 1;
    trials = {};
    for ii = 1:length(Spikes.Clusters)
        if ~isempty(Spikes.Clusters(ii).spikeTime)
            for i = 1:Behaviour.nCueMiss
                st = floor(Behaviour.cueMissTrace(i).LFPtime(1));
                stp = floor(Behaviour.cueMissTrace(i).LFPtime(end));
                temp = zeros(1,length(0:0.001:(stp-st)));
                spiketm = Spikes.Clusters(ii).spikeTime(Spikes.Clusters(ii).spikeTime>=st & Spikes.Clusters(ii).spikeTime<=stp);
                spiked = discretize(spiketm,st:0.001:stp); % 0.001 bin for 1 ms
                if ~isnan(spiked)
                    temp(spiked) = 1;
                end
                trials{count}(i,:) = temp;
            end
        else
            trials{count} = [];
        end
        count = count+1;
    end
    output = make_nice_mean_raster(trials,20,0);
    Spikes.PSTH.miss.spks = trials;
    Spikes.PSTH.miss.spkRates = output;
end

if isfield(Behaviour,'MIHitTrace')
    fprintf('Analyzing MI Hit trials...\n')
    count = 1;
    trials = {};
    for ii = 1:length(Spikes.Clusters)
        if ~isempty(Spikes.Clusters(ii).spikeTime)
            for i = 1:length(Behaviour.MIHitTrace)
                st = floor(Behaviour.MIHitTrace(i).LFPtime(1));
                stp = floor(Behaviour.MIHitTrace(i).LFPtime(end));
                temp = zeros(1,length(0:0.001:(stp-st)));
                spiketm = Spikes.Clusters(ii).spikeTime(Spikes.Clusters(ii).spikeTime>=st & Spikes.Clusters(ii).spikeTime<=stp);
                spiked = discretize(spiketm,st:0.001:stp); % 0.001 bin for 1 ms
                if ~isnan(spiked)
                    temp(spiked) = 1;
                end
                trials{count}(i,:) = temp;
            end
        else
            trials{count} = [];
        end
        count = count+1;
    end
    output = make_nice_mean_raster(trials,20,0);
    Spikes.PSTH.MIHit.spks = trials;
    Spikes.PSTH.MIHit.spkRates = output;
end

if isfield(Behaviour,'MIFATrace')
    fprintf('Analyzing MI false alarm trials...\n')
    count = 1;
    trials = {};
    for ii = 1:length(Spikes.Clusters)
        if ~isempty(Spikes.Clusters(ii).spikeTime)
            for i = 1:length(Behaviour.MIFATrace)
                st = floor(Behaviour.MIFATrace(i).LFPtime(1));
                stp = floor(Behaviour.MIFATrace(i).LFPtime(end));
                temp = zeros(1,length(0:0.001:(stp-st)));
                spiketm = Spikes.Clusters(ii).spikeTime(Spikes.Clusters(ii).spikeTime>=st & Spikes.Clusters(ii).spikeTime<=stp);
                spiked = discretize(spiketm,st:0.001:stp); % 0.001 bin for 1 ms
                if ~isnan(spiked)
                    temp(spiked) = 1;
                end
                trials{count}(i,:) = temp;
            end
        else
            trials{count} = [];
        end
        count = count+1;
    end
    output = make_nice_mean_raster(trials,20,0);
    Spikes.PSTH.MIFA.spks = trials;
    Spikes.PSTH.MIFA.spkRates = output;
end
%% Basic functions
    function output = make_nice_mean_raster(spmat,smooth_window,showplot)
        %*********** spmat1 and spmat2 are spike matrices of two conditions you wish to compare
        %*********** smooth_window ... gaussian smoothing in millisecs
        numconds = size(spmat,2);
        if (numconds==2)
            colo = [[1,0,0];[0,0,1]];
        else
            colo = jet(numconds);
        end
        for k = 1:numconds
            spud = spmat{k};
            numtrials = size(spud,1);
            smorate = gauss_smooth(sum( spud(1:numtrials,:))/....
                numtrials,smooth_window)*1000;
            if showplot
                plot(smorate,'k'); hold on;
                %                 set(H,'Color',colo(k,:));
            end
            output(k,:) = smorate;
            
        end
    end

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
    end
end