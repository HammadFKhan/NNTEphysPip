%% Load in avi data and calculate motion energy

motionEnergy = processVideoMotionEnergy();
%%

motionEnergyGUI()
%%
function motionEnergy = processVideoMotionEnergy()
% PROCESSVIDEOMOTIONENERGY - Load large .avi video in chunks and compute motion energy
% Motion energy = absolute difference between consecutive frames (as images)
% Uses incremental disk writing to avoid memory overflow
%
% Output:
%   motionEnergy - matfile object pointing to saved motion energy data

%% Step 1: User selects video file
[filename, filepath] = uigetfile('*.mp4', 'Select AVI Video File');

if isequal(filename, 0)
    disp('User canceled file selection');
    return;
end

fullFilePath = fullfile(filepath, filename);
fprintf('Selected file: %s\n', fullFilePath);

%% Step 2: Initialize video reader
vidObj = VideoReader(fullFilePath);

% Get video properties
totalFrames = vidObj.NumFrames;
frameHeight = vidObj.Height;
frameWidth = vidObj.Width;
frameRate = vidObj.FrameRate;

fprintf('Video properties:\n');
fprintf('  Total frames: %d\n', totalFrames);
fprintf('  Resolution: %d x %d\n', frameWidth, frameHeight);
fprintf('  Frame rate: %.2f fps\n', frameRate);
fprintf('  Approximate duration: %.2f seconds\n', totalFrames/frameRate);

%% Step 3: Setup output directory structure
baseDir = 'Y:\Hammad\Ephys\LeverTask\Data_for_Figures\Rebuttel\OrofacialData';
[~, videoName, ~] = fileparts(filename);

% Create video-specific folder
outputDir = fullfile(baseDir, videoName);
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created directory: %s\n', outputDir);
else
    fprintf('Using existing directory: %s\n', outputDir);
end

% Create output file path
saveName = fullfile(outputDir, sprintf('%s_motionEnergy.mat', videoName));

%% Step 4: Create output file and initialize with proper dimensions
% Create matfile object with write access
% Initialize the motion energy array in the file
motionEnergyTemp = zeros(frameHeight, frameWidth, totalFrames - 1, 'uint8');
save(saveName, 'motionEnergyTemp', 'frameRate', 'totalFrames', '-v7.3');
clear motionEnergyTemp; % Free memory immediately

% Now create matfile object for incremental writing
matObj = matfile(saveName, 'Writable', true);

fprintf('Created output file: %s\n', saveName);

%% Step 5: Set chunk size (adjust based on available memory)
chunkSize = 1000; % Process 500 frames at a time
numChunks = ceil(totalFrames / chunkSize);

%% Step 6: Process video in chunks
fprintf('\nProcessing video in %d chunks...\n', numChunks);

previousFrame = []; % Store last frame from previous chunk for continuity
meIdx = 1; % Motion energy frame index

for chunkIdx = 1:numChunks
    % Calculate frame range for this chunk
    startFrame = (chunkIdx - 1) * chunkSize + 1;
    endFrame = min(chunkIdx * chunkSize, totalFrames);
    numFramesInChunk = endFrame - startFrame + 1;

    fprintf('Processing chunk %d/%d (frames %d-%d)...\n', ...
        chunkIdx, numChunks, startFrame, endFrame);

    % Read chunk of frames
    vidObj.CurrentTime = (startFrame - 1) / frameRate;
    frames = zeros(frameHeight, frameWidth, numFramesInChunk, 'uint8');

    for i = 1:numFramesInChunk
        if hasFrame(vidObj)
            frame = readFrame(vidObj);
            % Convert to grayscale if RGB
            if size(frame, 3) == 3
                frame = rgb2gray(frame);
            end
            frames(:, :, i) = frame;
        end
    end

    %% Step 7: Compute motion energy for this chunk
    % Calculate how many motion energy frames we'll generate
    if chunkIdx > 1
        numMEFrames = numFramesInChunk; % Including boundary frame
    else
        numMEFrames = numFramesInChunk - 1; % First chunk has one less
    end

    % Preallocate temporary buffer for this chunk's motion energy
    meChunk = zeros(frameHeight, frameWidth, numMEFrames, 'uint8');
    chunkMEIdx = 1;

    % Handle boundary with previous chunk
    if chunkIdx > 1
        % Compute motion energy between last frame of previous chunk and first of current
        diffFrame = abs(int16(frames(:,:,1)) - int16(previousFrame));
        meChunk(:, :, chunkMEIdx) = uint8(diffFrame);
        chunkMEIdx = chunkMEIdx + 1;
        startIdx = 2;
    else
        startIdx = 2;
    end

    % Compute motion energy within chunk
    for i = startIdx:numFramesInChunk
        frame1 = frames(:, :, i-1);
        frame2 = frames(:, :, i);

        % Motion energy = absolute difference between frames
        diffFrame = abs(int16(frame2) - int16(frame1));
        meChunk(:, :, chunkMEIdx) = uint8(diffFrame);
        chunkMEIdx = chunkMEIdx + 1;
    end

    %% Step 8: Write chunk to disk
    % Write this chunk's motion energy to the matfile
    endMEIdx = meIdx + numMEFrames - 1;
    matObj.motionEnergyTemp(:, :, meIdx:endMEIdx) = meChunk;

    fprintf('  Written frames %d-%d to disk\n', meIdx, endMEIdx);

    % Update index for next chunk
    meIdx = endMEIdx + 1;

    % Store last frame for next chunk boundary
    previousFrame = frames(:, :, end);

    % Clear chunk data to free memory
    clear frames meChunk;

    % Progress update
    fprintf('  Completed %.1f%%\n', (chunkIdx/numChunks)*100);
end

fprintf('\nProcessing complete!\n');
fprintf('Motion energy dimensions: %d x %d x %d\n', ...
    frameHeight, frameWidth, totalFrames - 1);
fprintf('Results saved to: %s\n', saveName);

%% Step 9: Rename variable in file for cleaner output
% Load, rename, and resave
load(saveName, 'motionEnergyTemp');
motionEnergy = motionEnergyTemp;
clear motionEnergyTemp;
save(saveName, 'motionEnergy', 'frameRate', 'totalFrames', '-v7.3');

fprintf('Finalized motion energy array in file\n');

end


function motionEnergyGUI()
% MOTIONENERGYGUI - Interactive GUI for visualizing motion energy data
% 
% Fixed version with:
%   - ROI zoom updates during playback
%   - Proper cleanup of ROI lines when clearing

    %% Create main figure
    fig = uifigure('Name', 'Motion Energy Viewer', ...
                   'Position', [50 50 1600 900]);
    
    % Create grid layout (4 rows, 4 columns)
    mainGrid = uigridlayout(fig, [4 4]);
    mainGrid.RowHeight = {40, '1.5x', '1x', 35};
    mainGrid.ColumnWidth = {'1x', '1x', '1x', '1x'};
    
    %% Top panel - Control buttons
    controlPanel = uipanel(mainGrid);
    controlPanel.Layout.Row = 1;
    controlPanel.Layout.Column = [1 4];
    
    controlGrid = uigridlayout(controlPanel, [1 8]);
    controlGrid.ColumnWidth = {100, 100, 100, 120, 120, 120, '1x', 150};
    
    loadBtn = uibutton(controlGrid, 'Text', 'Load Data', ...
                       'ButtonPushedFcn', @loadDataCallback);
    loadBtn.Layout.Column = 1;
    
    playBtn = uibutton(controlGrid, 'Text', 'Play', ...
                       'ButtonPushedFcn', @playCallback, ...
                       'Enable', 'off');
    playBtn.Layout.Column = 2;
    
    stopBtn = uibutton(controlGrid, 'Text', 'Stop', ...
                       'ButtonPushedFcn', @stopCallback, ...
                       'Enable', 'off');
    stopBtn.Layout.Column = 3;
    
    addROIBtn = uibutton(controlGrid, 'Text', 'Add ROI', ...
                         'ButtonPushedFcn', @addROICallback, ...
                         'Enable', 'off');
    addROIBtn.Layout.Column = 4;
    
    computeROIBtn = uibutton(controlGrid, 'Text', 'Compute ROIs', ...
                             'ButtonPushedFcn', @computeROIsCallback, ...
                             'Enable', 'off');
    computeROIBtn.Layout.Column = 5;
    
    exportBtn = uibutton(controlGrid, 'Text', 'Export ROIs', ...
                         'ButtonPushedFcn', @exportROIsCallback, ...
                         'Enable', 'off');
    exportBtn.Layout.Column = 6;
    
    statusLabel = uilabel(controlGrid, 'Text', 'No data loaded');
    statusLabel.Layout.Column = 7;
    
    frameLabel = uilabel(controlGrid, 'Text', 'Frame: 0/0', ...
                        'HorizontalAlignment', 'right');
    frameLabel.Layout.Column = 8;
    
    %% Row 2, Col 1 - Motion energy image display
    imagePanel = uipanel(mainGrid, 'Title', 'Motion Energy Frame');
    imagePanel.Layout.Row = 2;
    imagePanel.Layout.Column = 1;
    
    imageGrid = uigridlayout(imagePanel, [2 1]);
    imageGrid.RowHeight = {'1x', 40};
    
    ax1 = uiaxes(imageGrid);
    ax1.Layout.Row = 1;
    axis(ax1, 'image');
    colormap(ax1, 'hot');
    title(ax1, 'Motion Energy');
    
    frameSlider = uislider(imageGrid, 'Limits', [1 100], ...
                          'Value', 1, ...
                          'ValueChangedFcn', @frameSliderCallback, ...
                          'Enable', 'off');
    frameSlider.Layout.Row = 2;
    
    %% Row 2, Col 2 - Mean motion projection
    projPanel = uipanel(mainGrid, 'Title', 'Mean Motion Projection');
    projPanel.Layout.Row = 2;
    projPanel.Layout.Column = 2;
    
    axProj = uiaxes(projPanel);
    axis(axProj, 'image');
    colormap(axProj, 'hot');
    title(axProj, 'Time-Averaged Motion Map');
    
    %% Row 2, Col 3 - Zoomed ROI view
    roiZoomPanel = uipanel(mainGrid, 'Title', 'Selected ROI');
    roiZoomPanel.Layout.Row = 2;
    roiZoomPanel.Layout.Column = 3;
    
    roiZoomGrid = uigridlayout(roiZoomPanel, [2 1]);
    roiZoomGrid.RowHeight = {'1x', 30};
    
    axZoom = uiaxes(roiZoomGrid);
    axZoom.Layout.Row = 1;
    axis(axZoom, 'image');
    colormap(axZoom, 'hot');
    title(axZoom, 'No ROI Selected');
    
    roiListBox = uilistbox(roiZoomGrid, ...
                          'Items', {}, ...
                          'ValueChangedFcn', @roiSelectionCallback);
    roiListBox.Layout.Row = 2;
    
    %% Row 2, Col 4 - ROI list and controls
    roiControlPanel = uipanel(mainGrid, 'Title', 'ROI Management');
    roiControlPanel.Layout.Row = 2;
    roiControlPanel.Layout.Column = 4;
    
    roiControlGrid = uigridlayout(roiControlPanel, [4 1]);
    roiControlGrid.RowHeight = {'1x', 40, 40, 40};
    
    roiTable = uitable(roiControlGrid, ...
                      'ColumnName', {'ROI', 'Color', 'Computed'}, ...
                      'ColumnWidth', {60, 80, 80}, ...
                      'RowName', {}, ...
                      'CellSelectionCallback', @roiTableSelectionCallback);
    roiTable.Layout.Row = 1;
    
    renameROIBtn = uibutton(roiControlGrid, 'Text', 'Rename ROI', ...
                           'ButtonPushedFcn', @renameROICallback, ...
                           'Enable', 'off');
    renameROIBtn.Layout.Row = 2;
    
    deleteROIBtn = uibutton(roiControlGrid, 'Text', 'Delete ROI', ...
                           'ButtonPushedFcn', @deleteROICallback, ...
                           'Enable', 'off');
    deleteROIBtn.Layout.Row = 3;
    
    clearROIBtn = uibutton(roiControlGrid, 'Text', 'Clear All ROIs', ...
                          'ButtonPushedFcn', @clearAllROIsCallback, ...
                          'Enable', 'off');
    clearROIBtn.Layout.Row = 4;
    
    %% Row 3 - Time series plots
    plotPanel1 = uipanel(mainGrid, 'Title', 'Average Motion Energy');
    plotPanel1.Layout.Row = 3;
    plotPanel1.Layout.Column = [1 2];
    
    ax2 = uiaxes(plotPanel1);
    title(ax2, 'Average Motion Energy Over Time');
    xlabel(ax2, 'Time (s)');
    ylabel(ax2, 'Mean Motion Energy');
    hold(ax2, 'on');
    
    plotPanel2 = uipanel(mainGrid, 'Title', 'ROI Motion Energy (10% Subset)');
    plotPanel2.Layout.Row = 3;
    plotPanel2.Layout.Column = [3 4];
    
    ax3 = uiaxes(plotPanel2);
    title(ax3, 'ROI Motion Energy');
    xlabel(ax3, 'Time (s)');
    ylabel(ax3, 'Mean Motion Energy');
    hold(ax3, 'on');
    
    %% Row 4 - Progress bar
    progressPanel = uipanel(mainGrid);
    progressPanel.Layout.Row = 4;
    progressPanel.Layout.Column = [1 4];
    
    progressGrid = uigridlayout(progressPanel, [1 2]);
    progressGrid.ColumnWidth = {'1x', 100};
    
    progressLabel = uilabel(progressGrid, 'Text', 'Ready');
    progressLabel.Layout.Column = 1;
    
    progressPercent = uilabel(progressGrid, 'Text', '', ...
                             'HorizontalAlignment', 'right');
    progressPercent.Layout.Column = 2;
    
    %% Data storage
    data = struct();
    data.matObj = [];
    data.filepath = '';
    data.filename = '';
    data.frameRate = [];
    data.totalFrames = 0;
    data.currentFrame = 1;
    data.isPlaying = false;
    data.rois = {};
    data.roiNames = {};
    data.roiColors = {};
    data.roiMotionEnergySubset = {};
    data.roiComputed = {};
    data.avgMotion = [];
    data.timeVector = [];
    data.numSubsetFrames = 0;
    data.selectedROI = [];
    data.selectedROIIdx = [];
    data.meanProjection = [];
    data.imageHandle = [];
    data.zoomImageHandle = []; % Handle for ROI zoom image
    data.roiTimeLines = [];
    
    %% Define ROI colors
    roiColorPalette = {
        [1 0.5 0],    % Orange
        [1 0 0],      % Red
        [1 1 0],      % Yellow
        [0 1 0],      % Green
        [0 1 1],      % Cyan
        [0 0 1],      % Blue
        [1 0 1],      % Magenta
        [0.5 0 0.5]   % Purple
    };
    
    %% Callback functions
    
    function loadDataCallback(~, ~)
        [filename, filepath] = uigetfile('*.mat', 'Select Motion Energy MAT File');
        
        if isequal(filename, 0)
            return;
        end
        
        fullPath = fullfile(filepath, filename);
        data.filepath = filepath;
        data.filename = filename;
        
        try
            progressLabel.Text = 'Opening file...';
            progressPercent.Text = '0%';
            drawnow;
            
            data.matObj = matfile(fullPath);
            data.frameRate = data.matObj.frameRate;
            data.totalFrames = data.matObj.totalFrames - 1;
            data.currentFrame = 1;
            
            data.numSubsetFrames = round(data.totalFrames * 0.1);
            
            progressLabel.Text = sprintf('Loading %d frames (10%%)...', data.numSubsetFrames);
            progressPercent.Text = '10%';
            drawnow;
            
            data.motionEnergySubset = data.matObj.motionEnergy(:, :, 1:data.numSubsetFrames);
            data.timeVectorSubset = (1:data.numSubsetFrames) / data.frameRate;
            
            progressLabel.Text = 'Computing mean projection...';
            progressPercent.Text = '50%';
            drawnow;
            
            data.meanProjection = mean(data.motionEnergySubset, 3);
            
            imagesc(axProj, data.meanProjection);
            colorbar(axProj);
            axis(axProj, 'image');
            title(axProj, 'Mean Motion Projection (First 10%)');
            
            progressLabel.Text = 'Computing average motion energy...';
            progressPercent.Text = '75%';
            drawnow;
            
            data.avgMotionSubset = squeeze(mean(mean(data.motionEnergySubset, 1), 2));
            
            playBtn.Enable = 'on';
            stopBtn.Enable = 'on';
            addROIBtn.Enable = 'on';
            frameSlider.Enable = 'on';
            frameSlider.Limits = [1 data.numSubsetFrames];
            frameSlider.Value = 1;
            
            frame = data.motionEnergySubset(:, :, 1);
            data.imageHandle = imagesc(ax1, frame);
            colorbar(ax1);
            axis(ax1, 'image');
            title(ax1, 'Motion Energy - Frame 1');
            
            cla(ax2);
            plot(ax2, data.timeVectorSubset, data.avgMotionSubset, 'k-', 'LineWidth', 1.5);
            data.timeLine = xline(ax2, data.timeVectorSubset(1), 'r', 'LineWidth', 2);
            
            windowSize = min(60, data.timeVectorSubset(end));
            xlim(ax2, [0 windowSize]);
            
            progressLabel.Text = sprintf('Ready | %s', filename);
            progressPercent.Text = '100%';
            statusLabel.Text = sprintf('Loaded: %s', filename);
            frameLabel.Text = sprintf('Frame: 1/%d (0.00s)', data.numSubsetFrames);
            
        catch ME
            uialert(fig, sprintf('Error loading file: %s', ME.message), 'Load Error');
            progressLabel.Text = 'Error loading data';
            progressPercent.Text = '';
        end
    end
    
    function displayFrame(frameNum, updateROIZoom)
        % updateROIZoom: whether to update the ROI zoom panel
        if isempty(data.motionEnergySubset)
            return;
        end
        
        frame = data.motionEnergySubset(:, :, frameNum);
        if ~isempty(data.imageHandle) && isvalid(data.imageHandle)
            data.imageHandle.CData = frame;
        else
            data.imageHandle = imagesc(ax1, frame);
            colorbar(ax1);
            axis(ax1, 'image');
        end
        
        title(ax1, sprintf('Motion Energy - Frame %d', frameNum));
        
        currentTime = data.timeVectorSubset(frameNum);
        if isfield(data, 'timeLine') && isvalid(data.timeLine)
            data.timeLine.Value = currentTime;
        end
        
        if ~isempty(data.roiTimeLines)
            for i = 1:length(data.roiTimeLines)
                if isvalid(data.roiTimeLines(i))
                    data.roiTimeLines(i).Value = currentTime;
                end
            end
        end
        
        windowSize = 60;
        if currentTime > windowSize/2
            xlim(ax2, [currentTime - windowSize/2, currentTime + windowSize/2]);
            xlim(ax3, [currentTime - windowSize/2, currentTime + windowSize/2]);
        else
            xlim(ax2, [0, windowSize]);
            xlim(ax3, [0, windowSize]);
        end
        
        frameLabel.Text = sprintf('Frame: %d/%d (%.2fs)', ...
                                  frameNum, data.numSubsetFrames, currentTime);
        
        if updateROIZoom
            drawROIsOnAxes();
            updateZoomedROI(frameNum);
        end
    end
    
    function drawROIsOnAxes()
        warnState = warning('off', 'all');
        
        hold(ax1, 'on');
        hold(axProj, 'on');
        
        for i = 1:length(data.rois)
            if isvalid(data.rois{i})
                try
                    data.rois{i}.Color = data.roiColors{i};
                    data.rois{i}.LineWidth = 2;
                    
                    pos = data.rois{i}.Position;
                    plot(axProj, pos(:,1), pos(:,2), ...
                         'Color', data.roiColors{i}, 'LineWidth', 2);
                catch
                    % Ignore drawing errors
                end
            end
        end
        
        hold(ax1, 'off');
        hold(axProj, 'off');
        
        warning(warnState);
    end
    
    function updateZoomedROI(frameNum)
        % Update ROI zoom panel with current frame
        if isempty(data.selectedROI) || data.selectedROI > length(data.rois)
            cla(axZoom);
            data.zoomImageHandle = [];
            title(axZoom, 'No ROI Selected');
            return;
        end
        
        roi = data.rois{data.selectedROI};
        if ~isvalid(roi)
            return;
        end
        
        frame = data.motionEnergySubset(:, :, frameNum);
        
        pos = roi.Position;
        xmin = max(1, floor(min(pos(:,1))) - 10);
        xmax = min(size(frame,2), ceil(max(pos(:,1))) + 10);
        ymin = max(1, floor(min(pos(:,2))) - 10);
        ymax = min(size(frame,1), ceil(max(pos(:,2))) + 10);
        
        % Extract zoomed region
        zoomedFrame = frame(ymin:ymax, xmin:xmax);
        
        % Update image efficiently (like main frame)
        if ~isempty(data.zoomImageHandle) && isvalid(data.zoomImageHandle)
            data.zoomImageHandle.CData = zoomedFrame;
        else
            data.zoomImageHandle = imagesc(axZoom, zoomedFrame);
            colorbar(axZoom);
            axis(axZoom, 'image');
        end
        
        title(axZoom, sprintf('%s Zoom', data.roiNames{data.selectedROI}));
        
        % Draw ROI outline
        hold(axZoom, 'on');
        adjustedPos = pos - [xmin ymin];
        plot(axZoom, adjustedPos(:,1), adjustedPos(:,2), ...
             'Color', data.roiColors{data.selectedROI}, 'LineWidth', 2);
        hold(axZoom, 'off');
    end
    
    function frameSliderCallback(src, ~)
        data.currentFrame = round(src.Value);
        displayFrame(data.currentFrame, true);
        drawnow;
    end
    
    function playCallback(~, ~)
        % NON-BLOCKING playback with ROI zoom updates
        data.isPlaying = true;
        playBtn.Enable = 'off';
        stopBtn.Enable = 'on';
        
        playbackSpeed = 1;
        targetFPS = data.frameRate * playbackSpeed;
        frameInterval = 1 / targetFPS;
        
        displayEveryNFrames = max(1, round(targetFPS / 30));
        
        frameCount = 0;
        lastTime = tic;
        
        while data.isPlaying && data.currentFrame <= data.numSubsetFrames
            % NON-BLOCKING wait
            if toc(lastTime) >= frameInterval
                
                if mod(frameCount, displayEveryNFrames) == 0
                    % Update both main frame AND ROI zoom during playback
                    displayFrame(data.currentFrame, true);
                    frameSlider.Value = data.currentFrame;
                end
                
                data.currentFrame = data.currentFrame + 1;
                frameCount = frameCount + 1;
                lastTime = tic;
            end
            
            % Process UI events (allows stop button to work)
            drawnow limitrate;
            
            if ~data.isPlaying
                break;
            end
        end
        
        data.isPlaying = false;
        playBtn.Enable = 'on';
        stopBtn.Enable = 'on';
        
        if data.currentFrame > data.numSubsetFrames
            data.currentFrame = 1;
            windowSize = min(60, data.timeVectorSubset(end));
            xlim(ax2, [0 windowSize]);
            xlim(ax3, [0 windowSize]);
        end
        
        % Final update
        displayFrame(data.currentFrame, true);
        drawnow;
    end
    
    function stopCallback(~, ~)
        data.isPlaying = false;
    end
    
    function addROICallback(~, ~)
        warnState = warning('off', 'all');
        
        choice = uiconfirm(fig, ...
            'Where would you like to draw the ROI?', ...
            'Select Drawing Canvas', ...
            'Options', {'Current Frame', 'Mean Projection', 'Cancel'}, ...
            'DefaultOption', 2);
        
        if strcmp(choice, 'Cancel')
            warning(warnState);
            return;
        end
        
        statusLabel.Text = 'Draw ROI...';
        
        if strcmp(choice, 'Mean Projection')
            drawAxis = axProj;
        else
            drawAxis = ax1;
        end
        
        roiIdx = length(data.rois) + 1;
        colorIdx = mod(roiIdx - 1, length(roiColorPalette)) + 1;
        roiColor = roiColorPalette{colorIdx};
        
        try
            roi = drawpolygon(drawAxis, 'Color', roiColor, 'LineWidth', 2);
            
            if ~isempty(roi)
                data.rois{end+1} = roi;
                data.roiNames{end+1} = sprintf('ROI%d', roiIdx);
                data.roiColors{end+1} = roiColor;
                data.roiMotionEnergySubset{end+1} = [];
                data.roiComputed{end+1} = false;
                
                roiListBox.Items = data.roiNames;
                updateROITable();
                
                computeROIBtn.Enable = 'on';
                renameROIBtn.Enable = 'on';
                deleteROIBtn.Enable = 'on';
                clearROIBtn.Enable = 'on';
                
                displayFrame(data.currentFrame, true);
                drawnow;
                
                statusLabel.Text = sprintf('Added %s', data.roiNames{end});
            end
        catch ME
            uialert(fig, sprintf('Error drawing ROI: %s', ME.message), 'ROI Error');
        end
        
        warning(warnState);
    end
    
    function renameROICallback(~, ~)
        if isempty(data.selectedROIIdx)
            uialert(fig, 'Please select an ROI from the table first.', 'No Selection');
            return;
        end
        
        idx = data.selectedROIIdx;
        currentName = data.roiNames{idx};
        
        answer = inputdlg({'Enter new ROI name:'}, 'Rename ROI', 1, {currentName});
        
        if ~isempty(answer)
            data.roiNames{idx} = answer{1};
            roiListBox.Items = data.roiNames;
            updateROITable();
            statusLabel.Text = sprintf('Renamed to %s', answer{1});
        end
    end
    
    function computeROIsCallback(~, ~)
        if isempty(data.rois)
            return;
        end
        
        cla(ax3);
        hold(ax3, 'on');
        
        data.roiTimeLines = [];
        
        for roiIdx = 1:length(data.rois)
            roi = data.rois{roiIdx};
            if ~isvalid(roi)
                continue;
            end
            
            progressLabel.Text = sprintf('Computing %s on subset (%d/%d)...', ...
                                         data.roiNames{roiIdx}, roiIdx, length(data.rois));
            progressPercent.Text = sprintf('%.0f%%', (roiIdx-1)/length(data.rois)*100);
            drawnow;
            
            firstFrame = data.motionEnergySubset(:,:,1);
            [h, w] = size(firstFrame);
            roiMask = poly2mask(roi.Position(:,1), roi.Position(:,2), h, w);
            
            roiMotionSubset = zeros(data.numSubsetFrames, 1);
            for i = 1:data.numSubsetFrames
                frame = double(data.motionEnergySubset(:, :, i));
                roiMotionSubset(i) = mean(frame(roiMask));
            end
            
            data.roiMotionEnergySubset{roiIdx} = roiMotionSubset;
            data.roiComputed{roiIdx} = true;
            
            plot(ax3, data.timeVectorSubset, roiMotionSubset, ...
                 'Color', data.roiColors{roiIdx}, 'LineWidth', 1.5, ...
                 'DisplayName', data.roiNames{roiIdx});
        end
        
        if ~isempty(data.roiMotionEnergySubset)
            currentTime = data.timeVectorSubset(data.currentFrame);
            data.roiTimeLines = xline(ax3, currentTime, 'r', 'LineWidth', 2);
            
            windowSize = min(60, data.timeVectorSubset(end));
            xlim(ax3, [0, windowSize]);
        end
        
        legend(ax3, 'Location', 'best');
        hold(ax3, 'off');
        
        updateROITable();
        exportBtn.Enable = 'on';
        
        progressLabel.Text = 'ROI computation complete (subset only - export for full data)';
        progressPercent.Text = '100%';
        statusLabel.Text = 'ROI computed on 10% subset';
    end
    
    function exportROIsCallback(~, ~)
        if isempty(data.rois)
            uialert(fig, 'No ROIs to export.', 'No ROIs');
            return;
        end
        
        defaultName = strrep(data.filename, '_motionEnergy.mat', '_ROI_data.mat');
        [outFilename, outFilepath] = uiputfile('*.mat', 'Save ROI Motion Energy', ...
                                               fullfile(data.filepath, defaultName));
        
        if isequal(outFilename, 0)
            return;
        end
        
        outputPath = fullfile(outFilepath, outFilename);
        
        try
            progressLabel.Text = 'Computing full dataset ROIs...';
            drawnow;
            
            exportData = struct();
            exportData.frameRate = data.frameRate;
            exportData.totalFrames = data.totalFrames;
            exportData.timeVector = (1:data.totalFrames)' / data.frameRate;
            
            for i = 1:length(data.roiNames)
                exportData.([data.roiNames{i} '_Position']) = data.rois{i}.Position;
            end
            
            save(outputPath, '-struct', 'exportData', '-v7.3');
            matObj = matfile(outputPath, 'Writable', true);
            
            for roiIdx = 1:length(data.rois)
                roi = data.rois{roiIdx};
                if ~isvalid(roi)
                    continue;
                end
                
                roiName = data.roiNames{roiIdx};
                progressLabel.Text = sprintf('Computing %s on full dataset (%d/%d)...', ...
                                             roiName, roiIdx, length(data.rois));
                progressPercent.Text = sprintf('%.0f%%', (roiIdx-1)/length(data.rois)*100);
                drawnow;
                
                firstFrame = data.motionEnergySubset(:,:,1);
                [h, w] = size(firstFrame);
                roiMask = poly2mask(roi.Position(:,1), roi.Position(:,2), h, w);
                
                roiMotionFull = zeros(data.totalFrames, 1);
                
                chunkSize = 1000;
                numChunks = ceil(data.totalFrames / chunkSize);
                
                for chunkIdx = 1:numChunks
                    startIdx = (chunkIdx-1)*chunkSize + 1;
                    endIdx = min(chunkIdx*chunkSize, data.totalFrames);
                    
                    meChunk = data.matObj.motionEnergy(:, :, startIdx:endIdx);
                    
                    for j = 1:(endIdx - startIdx + 1)
                        frame = double(meChunk(:, :, j));
                        roiMotionFull(startIdx + j - 1) = mean(frame(roiMask));
                    end
                    
                    chunkProgress = chunkIdx / numChunks;
                    totalProgress = ((roiIdx-1) + chunkProgress) / length(data.rois);
                    progressPercent.Text = sprintf('%.0f%%', totalProgress*100);
                    
                    if mod(chunkIdx, 10) == 0
                        drawnow;
                    end
                end
                
                matObj.(roiName) = roiMotionFull;
                
                progressLabel.Text = sprintf('Saved %s (%d/%d)', ...
                                             roiName, roiIdx, length(data.rois));
                drawnow;
            end
            
            progressLabel.Text = 'Export complete!';
            progressPercent.Text = '100%';
            
            uialert(fig, sprintf('ROI data (full dataset) exported to:\n%s\n\nAll %d ROIs computed on full %d frames.', ...
                                 outputPath, length(data.rois), data.totalFrames), ...
                    'Export Successful', 'Icon', 'success');
            
            statusLabel.Text = sprintf('Exported to %s', outFilename);
            
        catch ME
            uialert(fig, sprintf('Error during export: %s', ME.message), 'Export Error');
            progressLabel.Text = 'Export failed';
        end
    end
    
    function deleteROICallback(~, ~)
        if isempty(data.selectedROIIdx)
            uialert(fig, 'Please select an ROI from the table first.', 'No Selection');
            return;
        end
        
        idx = data.selectedROIIdx;
        
        if isvalid(data.rois{idx})
            delete(data.rois{idx});
        end
        
        data.rois(idx) = [];
        data.roiNames(idx) = [];
        data.roiColors(idx) = [];
        data.roiMotionEnergySubset(idx) = [];
        data.roiComputed(idx) = [];
        
        roiListBox.Items = data.roiNames;
        updateROITable();
        
        data.selectedROI = [];
        data.selectedROIIdx = [];
        cla(axZoom);
        data.zoomImageHandle = [];
        title(axZoom, 'No ROI Selected');
        
        displayFrame(data.currentFrame, true);
        drawnow;
    end
    
    function clearAllROIsCallback(~, ~)
        % Delete all ROI objects
        for i = 1:length(data.rois)
            if isvalid(data.rois{i})
                delete(data.rois{i});
            end
        end
        
        % Clear all data
        data.rois = {};
        data.roiNames = {};
        data.roiColors = {};
        data.roiMotionEnergySubset = {};
        data.roiComputed = {};
        data.selectedROI = [];
        data.selectedROIIdx = [];
        
        % Clear UI elements
        roiListBox.Items = {};
        updateROITable();
        
        % Clear all axes properly
        cla(axZoom);
        data.zoomImageHandle = [];
        title(axZoom, 'No ROI Selected');
        
        % Clear ROI plot completely
        cla(ax3);
        hold(ax3, 'on');
        title(ax3, 'ROI Motion Energy');
        xlabel(ax3, 'Time (s)');
        ylabel(ax3, 'Mean Motion Energy');
        
        % Clear time lines
        data.roiTimeLines = [];
        
        % Disable buttons
        computeROIBtn.Enable = 'off';
        renameROIBtn.Enable = 'off';
        deleteROIBtn.Enable = 'off';
        clearROIBtn.Enable = 'off';
        exportBtn.Enable = 'off';
        
        % Redraw main frame
        displayFrame(data.currentFrame, true);
        drawnow;
    end
    
    function roiSelectionCallback(src, ~)
        selectedName = src.Value;
        if isempty(selectedName)
            return;
        end
        
        data.selectedROI = find(strcmp(data.roiNames, selectedName), 1);
        data.selectedROIIdx = data.selectedROI;
        displayFrame(data.currentFrame, true);
        drawnow;
    end
    
    function roiTableSelectionCallback(~, event)
        if ~isempty(event.Indices)
            data.selectedROIIdx = event.Indices(1);
            data.selectedROI = data.selectedROIIdx;
            
            if data.selectedROIIdx <= length(data.roiNames)
                roiListBox.Value = data.roiNames{data.selectedROIIdx};
                displayFrame(data.currentFrame, true);
                drawnow;
            end
        end
    end
    
    function updateROITable()
        if isempty(data.rois)
            roiTable.Data = {};
            return;
        end
        
        tableData = cell(length(data.rois), 3);
        for i = 1:length(data.rois)
            tableData{i, 1} = data.roiNames{i};
            tableData{i, 2} = sprintf('[%.1f %.1f %.1f]', data.roiColors{i});
            if data.roiComputed{i}
                tableData{i, 3} = '✓';
            else
                tableData{i, 3} = '';
            end
        end
        roiTable.Data = tableData;
    end

end
