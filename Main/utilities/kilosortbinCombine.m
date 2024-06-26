function fpath = kilosortbinCombine()
%% Function to combine kilsort bin files for multiple session spike sorting.
% All you need is a folder with two kilosort bin files (name starting with
% 'kilosort' (example. kilosortControl, kilosortOpto). Then select the bin
% files when prompted.

%Caution: Check Issue #132 on Mouseland kilosort for potential
%complications with this approach. TLDR: The correlation matrix during the
%spike sorting can generate a very large matrix if the drift correction
%fails. Causing memory bottlenecking on the GPU.

[fname,fpath] = uigetfile('*.bin','MultiSelect','on');
if length(fname)<2
    error('Only 1 file was provided')
end
fidW        = fopen(strcat(fpath,'kilosort.bin'),   'w+'); % open for writing processed data
for numFile = 1:length(fname)
    fbinary = fullfile(fpath, fname{numFile});
    trange    = [0 Inf]; % time range to sort
    bytes       = get_file_size(fbinary); % size in bytes of raw binary
    NchanTOT = 64;
    nTimepoints = floor(bytes/NchanTOT/2); % number of total timepoints
    tstart  = ceil(trange(1) * 20000); % starting timepoint for processing data segment
    tend    = min(nTimepoints, ceil(trange(2) * 20000)); % ending timepoint
    sampsToRead = tend-tstart; % total number of samples to read
    twind = tstart * NchanTOT*2; % skip this many bytes at the start
    NT = 1000000;
    Nbatch      = ceil(sampsToRead /NT); % number of data batches
    
    fid         = fopen(fbinary, 'r'); % open for reading raw data
    
    for ibatch = 1:Nbatch
        disp(['Batch: ' num2str(ibatch)])
        offset = 0;
        fseek(fid, offset, 'bof'); % fseek to batch start in raw file
        buff = fread(fid, [NchanTOT NT], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
        count = fwrite(fidW, buff, 'int16'); % write this batch to binary file
    end
    fclose(fid);
end
fclose(fidW);
% delete uneccesary bin files for kilosort to find
for numFile = 1:length(fname)
    fbinary = fullfile(fpath, fname{numFile});
    delete(fbinary)
end
disp('Kilosort file combined')
