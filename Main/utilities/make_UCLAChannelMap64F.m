function make_UCLAChannelMap64F(fpath)
% create a channel Map file for simulated data (eMouse)

% here I know a priori what order my channels are in.  So I just manually 
% make a list of channel indices (and give
% an index to dead channels too). chanMap(1) is the row in the raw binary
% file for the first channel. chanMap(1:2) = [33 34] in my case, which happen to
% be dead channels. 
load('UCLA_chanmap_64F2.mat')
% chanMap = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 ...
%     25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 ...
%     50 51 52 53 54 55 56 57 58 59 60 61 62 63 64];
chanMap = [1:64];
% the first thing Kilosort does is reorder the data with data = data(chanMap, :).
% Now we declare which channels are "connected" in this normal ordering, 
% meaning not dead or used for non-ephys data

connected = true(64, 1); 
% connected(33) = 0; %For no connections

% now we define the horizontal (x) and vertical (y) coordinates of these
% channels. For dead or nonephys channels the values won't matter. Again
% I will take this information from the specifications of the probe. These
% are in um here, but the absolute scaling doesn't really matter in the
% algorithm. 

% xcoords = ones(1,64);

xcoords = s.sorted_probe_wiring(:,2)';%[1*ones(1, 21), 2*ones(1, 22), 3*ones(1, 21)];
% ycoords = [15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 16 17 18 19 20 21 ...
%     22 20 18 16 14 12 10 8 6 4 2 1 3 5 7 9 11 13 15 17 19 21 ...
%     21 20 19 18 17 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
ycoords = s.sorted_probe_wiring(:,4)';

% Often, multi-shank probes or tetrodes will be organized into groups of
% channels that cannot possibly share spikes with the rest of the probe. This helps
% the algorithm discard noisy templates shared across groups. In
% this case, we set kcoords to indicate which group the channel belongs to.
% In our case all channels are on the same shank in a single group so we
% assign them all to group 1. 

kcoords = s.sorted_probe_wiring(:,5)';

% at this point in Kilosort we do data = data(connected, :), ycoords =
% ycoords(connected), xcoords = xcoords(connected) and kcoords =
% kcoords(connected) and no more channel map information is needed (in particular
% no "adjacency graphs" like in KlustaKwik). 
% Now we can save our channel map for the eMouse. 

% would be good to also save the sampling frequency here
fs = 20000; 

save(fullfile(fpath, 'chanMap64F2.mat'), 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs')