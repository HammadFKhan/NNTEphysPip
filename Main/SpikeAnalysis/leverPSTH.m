function Spikes = leverPSTH(Spikes,Behaviour)
% Organizes and calculates smoothed firing rates across trials
% This data is organized for coherence measurements
if ~isstruct(Behaviour),error('Behaviour must be structure!');end
if ~isstruct(Spikes), error('Spikes must be a structure!');end

%% Analyze for hit trials
if isfield(Behaviour,'nHit')
    fprintf('Analyzing hit trials...\n')
    count = 1;
    trials = {};
    for ii = 1:length(Spikes.Clusters)
        if ~isempty(Spikes.Clusters(ii).spikeTime)
            for i = 1:Behaviour.nHit
                temp = zeros(1,length(0:0.001:Behaviour.hitTrace(i).t2-Behaviour.hitTrace(i).t1));
                spiketm = Spikes.Clusters(ii).spikeTime(Spikes.Clusters(ii).spikeTime>=Behaviour.hitTrace(i).t1 & Spikes.Clusters(ii).spikeTime<=Behaviour.hitTrace(i).t2);
                spiked = discretize(spiketm,Behaviour.hitTrace(i).t1:0.001:Behaviour.hitTrace(i).t2); % 0.001 bin for 1 ms
                if ~isnan(spiked)
                    temp(spiked) = 1;
                end
                trials{count}(i,:) = temp;
                
            end
        else
            continue
        end
        count = count+1;
    end
    %%
    output = make_nice_mean_raster(trials,20,0);
    %%
    Spikes.PSTH.hit = trials;
    Spikes.PSTH.hitspikeRates = output;
end
%% Analyze for miss trials
if isfield(Behaviour,'nMiss')
    fprintf('Analyzing miss trials...\n')
    count = 1;
    trials = {};
    for ii = 1:length(Spikes.Clusters)
        if ~isempty(Spikes.Clusters(ii).spikeTime)
            for i = 1:Behaviour.nMiss
                temp= zeros(1,length(0:0.001:Behaviour.missTrace(i).t2-Behaviour.missTrace(i).t1));
                spiketm = Spikes.Clusters(ii).spikeTime(Spikes.Clusters(ii).spikeTime>=Behaviour.missTrace(i).t1 & Spikes.Clusters(ii).spikeTime<=Behaviour.missTrace(i).t2);
                spiked = discretize(spiketm,Behaviour.missTrace(i).t1:0.001:Behaviour.missTrace(i).t2);
                if ~isnan(spiked)
                    temp(spiked) = 1;
                end
                trials{count}(i,:) = temp;
            end
        else
            continue
        end
        count = count+1;
    end
    output = make_nice_mean_raster(trials,20,0);
    Spikes.PSTH.miss = trials;
    Spikes.PSTH.missspikeRates = output;
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