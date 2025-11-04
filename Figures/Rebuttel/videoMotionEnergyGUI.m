function videoMotionEnergyGUI()
% VIDEOMOTIONENERGYGUI - Body part segmentation with SVD motion energy
%
% Features:
%   - Predefined rectangular ROIs: Tongue, Limb, Whiskers, Body, Pupil
%   - SVD-based motion energy computation (whole-frame + per-ROI)
%   - Space-efficient output: motSVD + uMotMask + per-ROI SVD

%% Create main figure
fig = uifigure('Name', 'Video Motion Energy Processor', ...
    'Position', [50 50 1800 900]);

% Create grid layout (5 rows, 5 columns)
mainGrid = uigridlayout(fig, [5 5]);
mainGrid.RowHeight = {40, 40, '1.5x', '1x', 35};
mainGrid.ColumnWidth = {'1x', '1x', '1x', '1x', '1x'};

%% Row 1 - File management controls
filePanel = uipanel(mainGrid);
filePanel.Layout.Row = 1;
filePanel.Layout.Column = [1 5];

fileGrid = uigridlayout(filePanel, [1 5]);
fileGrid.ColumnWidth = {150, 150, 150, '1x', 200};

loadVideoBtn = uibutton(fileGrid, 'Text', 'Load Video', ...
    'ButtonPushedFcn', @loadVideoCallback);
loadVideoBtn.Layout.Column = 1;

setOutputDirBtn = uibutton(fileGrid, 'Text', 'Set Output Dir', ...
    'ButtonPushedFcn', @setOutputDirCallback);
setOutputDirBtn.Layout.Column = 2;

processBtn = uibutton(fileGrid, 'Text', 'Process Video', ...
    'ButtonPushedFcn', @processVideoCallback, ...
    'Enable', 'off');
processBtn.Layout.Column = 3;

outputDirLabel = uilabel(fileGrid, 'Text', 'Output: Not set');
outputDirLabel.Layout.Column = 4;

videoInfoLabel = uilabel(fileGrid, 'Text', 'No video loaded', ...
    'HorizontalAlignment', 'right');
videoInfoLabel.Layout.Column = 5;

%% Row 2 - Playback and ROI controls
controlPanel = uipanel(mainGrid);
controlPanel.Layout.Row = 2;
controlPanel.Layout.Column = [1 5];

controlGrid = uigridlayout(controlPanel, [1 11]);
controlGrid.ColumnWidth = {80, 80, 100, 100, 100, 100, 100, 120, 100, '1x', 150};

playBtn = uibutton(controlGrid, 'Text', 'Play', ...
    'ButtonPushedFcn', @playCallback, ...
    'Enable', 'off');
playBtn.Layout.Column = 1;

stopBtn = uibutton(controlGrid, 'Text', 'Stop', ...
    'ButtonPushedFcn', @stopCallback, ...
    'Enable', 'off');
stopBtn.Layout.Column = 2;

tongueROIBtn = uibutton(controlGrid, 'Text', 'Add Tongue', ...
    'ButtonPushedFcn', @(~,~) addBodyPartROI('Tongue'), ...
    'Enable', 'off', ...
    'BackgroundColor', [1 0.5 0]);
tongueROIBtn.Layout.Column = 3;

limbROIBtn = uibutton(controlGrid, 'Text', 'Add Limb', ...
    'ButtonPushedFcn', @(~,~) addBodyPartROI('Limb'), ...
    'Enable', 'off', ...
    'BackgroundColor', [1 0 0]);
limbROIBtn.Layout.Column = 4;

whiskersROIBtn = uibutton(controlGrid, 'Text', 'Add Whiskers', ...
    'ButtonPushedFcn', @(~,~) addBodyPartROI('Whiskers'), ...
    'Enable', 'off', ...
    'BackgroundColor', [0 1 0]);
whiskersROIBtn.Layout.Column = 5;

bodyROIBtn = uibutton(controlGrid, 'Text', 'Add Body', ...
    'ButtonPushedFcn', @(~,~) addBodyPartROI('Body'), ...
    'Enable', 'off', ...
    'BackgroundColor', [0 0 1]);
bodyROIBtn.Layout.Column = 6;

pupilROIBtn = uibutton(controlGrid, 'Text', 'Add Pupil', ...
    'ButtonPushedFcn', @(~,~) addBodyPartROI('Pupil'), ...
    'Enable', 'off', ...
    'BackgroundColor', [0 1 1]);
pupilROIBtn.Layout.Column = 7;

deleteROIBtn = uibutton(controlGrid, 'Text', 'Delete ROI', ...
    'ButtonPushedFcn', @deleteROICallback, ...
    'Enable', 'off');
deleteROIBtn.Layout.Column = 8;

clearROIBtn = uibutton(controlGrid, 'Text', 'Clear All', ...
    'ButtonPushedFcn', @clearAllROIsCallback, ...
    'Enable', 'off');
clearROIBtn.Layout.Column = 9;

statusLabel = uilabel(controlGrid, 'Text', 'Load video to begin');
statusLabel.Layout.Column = 10;

frameLabel = uilabel(controlGrid, 'Text', 'Frame: 0/0', ...
    'HorizontalAlignment', 'right');
frameLabel.Layout.Column = 11;

%% Row 3 - Video displays
rawVideoPanel = uipanel(mainGrid, 'Title', 'Raw Video Frame');
rawVideoPanel.Layout.Row = 3;
rawVideoPanel.Layout.Column = 1;

rawGrid = uigridlayout(rawVideoPanel, [2 1]);
rawGrid.RowHeight = {'1x', 40};

axRaw = uiaxes(rawGrid);
axRaw.Layout.Row = 1;
axis(axRaw, 'image');
colormap(axRaw, 'gray');
title(axRaw, 'Raw Video');

frameSlider = uislider(rawGrid, 'Limits', [1 100], ...
    'Value', 1, ...
    'ValueChangedFcn', @frameSliderCallback, ...
    'Enable', 'off');
frameSlider.Layout.Row = 2;

mePanel = uipanel(mainGrid, 'Title', 'Motion Energy Frame');
mePanel.Layout.Row = 3;
mePanel.Layout.Column = 2;

axME = uiaxes(mePanel);
axis(axME, 'image');
colormap(axME, 'hot');
title(axME, 'Motion Energy');

projPanel = uipanel(mainGrid, 'Title', 'Mean Motion Projection');
projPanel.Layout.Row = 3;
projPanel.Layout.Column = 3;

axProj = uiaxes(projPanel);
axis(axProj, 'image');
colormap(axProj, 'hot');
title(axProj, 'Time-Averaged Motion Map');

roiZoomPanel = uipanel(mainGrid, 'Title', 'Selected ROI (Raw)');
roiZoomPanel.Layout.Row = 3;
roiZoomPanel.Layout.Column = 4;

axZoom = uiaxes(roiZoomPanel);
axis(axZoom, 'image');
colormap(axZoom, 'gray');
title(axZoom, 'No ROI Selected');

roiControlPanel = uipanel(mainGrid, 'Title', 'Body Part ROIs');
roiControlPanel.Layout.Row = 3;
roiControlPanel.Layout.Column = 5;

roiControlGrid = uigridlayout(roiControlPanel, [2 1]);
roiControlGrid.RowHeight = {'1x', 30};

roiTable = uitable(roiControlGrid, ...
    'ColumnName', {'Body Part', 'Color'}, ...
    'ColumnWidth', {100, 80}, ...
    'RowName', {}, ...
    'CellSelectionCallback', @roiTableSelectionCallback);
roiTable.Layout.Row = 1;

roiListBox = uilistbox(roiControlGrid, ...
    'Items', {}, ...
    'ValueChangedFcn', @roiSelectionCallback);
roiListBox.Layout.Row = 2;

%% Row 4 - Time series plots
plotPanel1 = uipanel(mainGrid, 'Title', 'Average Motion Energy');
plotPanel1.Layout.Row = 4;
plotPanel1.Layout.Column = [1 3];

ax2 = uiaxes(plotPanel1);
title(ax2, 'Average Motion Energy Over Time');
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'Mean Motion Energy');
hold(ax2, 'on');

plotPanel2 = uipanel(mainGrid, 'Title', 'Body Part Motion Energy (5% Subset)');
plotPanel2.Layout.Row = 4;
plotPanel2.Layout.Column = [4 5];

ax3 = uiaxes(plotPanel2);
title(ax3, 'Body Part Motion Energy');
xlabel(ax3, 'Time (s)');
ylabel(ax3, 'Mean Motion Energy');
hold(ax3, 'on');

%% Row 5 - Progress bar
progressPanel = uipanel(mainGrid);
progressPanel.Layout.Row = 5;
progressPanel.Layout.Column = [1 5];

progressGrid = uigridlayout(progressPanel, [1 2]);
progressGrid.ColumnWidth = {'1x', 100};

progressLabel = uilabel(progressGrid, 'Text', 'Ready');
progressLabel.Layout.Column = 1;

progressPercent = uilabel(progressGrid, 'Text', '', ...
    'HorizontalAlignment', 'right');
progressPercent.Layout.Column = 2;

%% Data storage
data = struct();
data.videoPath = '';
data.videoName = '';
data.outputDir = '';
data.videoSubset = [];
data.meSubset = [];
data.meanProjection = [];
data.frameRate = [];
data.totalFrames = 0;
data.subsetFrames = 0;
data.currentFrame = 1;
data.isPlaying = false;
data.rois = {};
data.roiNames = {};
data.roiColors = {};
data.roiRects = {};
data.roiMotionEnergySubset = {};
data.selectedROI = [];
data.selectedROIIdx = [];
data.rawImageHandle = [];
data.meImageHandle = [];
data.zoomImageHandle = [];
data.timeLine = [];
data.roiTimeLines = [];

%% Body part colors
bodyPartColors = struct('Tongue', [1 0.5 0], ...
    'Limb', [1 0 0], ...
    'Whiskers', [0 1 0], ...
    'Body', [0 0 1], ...
    'Pupil', [0 1 1]);

%% Callback functions

    function loadVideoCallback(~, ~)
        [filename, filepath] = uigetfile('*.avi', 'Select Video File');

        if isequal(filename, 0)
            return;
        end

        data.videoPath = fullfile(filepath, filename);
        data.videoName = filename;

        loadVideo();
    end

    function loadVideo()
        try
            progressLabel.Text = sprintf('Loading %s...', data.videoName);
            progressPercent.Text = '0%';
            drawnow;

            vidObj = VideoReader(data.videoPath);
            data.frameRate = vidObj.FrameRate;
            data.totalFrames = vidObj.NumFrames;
            data.subsetFrames = round(data.totalFrames * 0.05);

            videoInfoLabel.Text = sprintf('%d frames @ %.1f fps', data.totalFrames, data.frameRate);

            progressLabel.Text = sprintf('Loading %d frames (5%%)...', data.subsetFrames);
            drawnow;

            frameHeight = vidObj.Height;
            frameWidth = vidObj.Width;
            data.videoSubset = zeros(frameHeight, frameWidth, data.subsetFrames, 'uint8');

            for i = 1:data.subsetFrames
                if hasFrame(vidObj)
                    frame = readFrame(vidObj);
                    if size(frame, 3) == 3
                        frame = rgb2gray(frame);
                    end
                    data.videoSubset(:, :, i) = frame;
                end

                if mod(i, 100) == 0
                    progressPercent.Text = sprintf('%.0f%%', (i/data.subsetFrames)*50);
                    drawnow limitrate;
                end
            end

            progressLabel.Text = 'Computing motion energy...';
            drawnow;

            data.meSubset = zeros(frameHeight, frameWidth, data.subsetFrames-1, 'uint8');
            for i = 2:data.subsetFrames
                frame1 = int16(data.videoSubset(:, :, i-1));
                frame2 = int16(data.videoSubset(:, :, i));
                data.meSubset(:, :, i-1) = uint8(abs(frame2 - frame1));
            end

            data.meanProjection = mean(data.meSubset, 3);
            imagesc(axProj, data.meanProjection);
            colorbar(axProj);
            axis(axProj, 'image');
            title(axProj, 'Mean Motion Projection (First 5%)');

            % Enable controls
            playBtn.Enable = 'on';
            stopBtn.Enable = 'on';
            tongueROIBtn.Enable = 'on';
            limbROIBtn.Enable = 'on';
            whiskersROIBtn.Enable = 'on';
            bodyROIBtn.Enable = 'on';
            pupilROIBtn.Enable = 'on';
            frameSlider.Enable = 'on';
            frameSlider.Limits = [1 data.subsetFrames];
            frameSlider.Value = 1;

            data.currentFrame = 1;

            displayFrame(1, true);

            % Plot average motion energy
            avgMotion = squeeze(mean(mean(data.meSubset, 1), 2));
            timeVector = (1:length(avgMotion)) / data.frameRate;

            cla(ax2);
            plot(ax2, timeVector, avgMotion, 'k-', 'LineWidth', 1.5);
            data.timeLine = xline(ax2, timeVector(1), 'r', 'LineWidth', 2);
            windowSize = min(60, timeVector(end));
            xlim(ax2, [0 windowSize]);

            processBtn.Enable = 'on';

            progressLabel.Text = sprintf('Ready | %s', data.videoName);
            progressPercent.Text = '100%';
            statusLabel.Text = sprintf('Loaded: %s', data.videoName);
            frameLabel.Text = sprintf('Frame: 1/%d', data.subsetFrames);

        catch ME
            uialert(fig, sprintf('Error loading video: %s', ME.message), 'Load Error');
            progressLabel.Text = 'Error loading video';
        end
    end

    function setOutputDirCallback(~, ~)
        selectedDir = uigetdir('', 'Select Output Directory');

        if isequal(selectedDir, 0)
            return;
        end

        data.outputDir = selectedDir;
        outputDirLabel.Text = sprintf('Output: %s', selectedDir);
        statusLabel.Text = 'Output directory set';
    end

    function addBodyPartROI(bodyPart)
        warnState = warning('off', 'all');

        try
            roi = drawrectangle(axProj, 'Color', bodyPartColors.(bodyPart), 'LineWidth', 2);

            if ~isempty(roi)
                data.rois{end+1} = roi;
                data.roiNames{end+1} = bodyPart;
                data.roiColors{end+1} = bodyPartColors.(bodyPart);
                data.roiRects{end+1} = roi.Position;
                data.roiMotionEnergySubset{end+1} = [];

                roiListBox.Items = data.roiNames;
                updateROITable();

                deleteROIBtn.Enable = 'on';
                clearROIBtn.Enable = 'on';

                computeROISubset(length(data.rois));

                displayFrame(data.currentFrame, true);
                drawnow;

                statusLabel.Text = sprintf('Added %s ROI', bodyPart);
            end
        catch ME
            uialert(fig, sprintf('Error drawing ROI: %s', ME.message), 'ROI Error');
        end

        warning(warnState);
    end

    function computeROISubset(roiIdx)
        roi = data.rois{roiIdx};
        rect = roi.Position;

        [h, w] = size(data.meSubset(:,:,1));
        roiMask = false(h, w);
        xStart = max(1, round(rect(1)));
        xEnd = min(w, round(rect(1) + rect(3)));
        yStart = max(1, round(rect(2)));
        yEnd = min(h, round(rect(2) + rect(4)));
        roiMask(yStart:yEnd, xStart:xEnd) = true;

        roiMotion = zeros(size(data.meSubset, 3), 1);
        for i = 1:size(data.meSubset, 3)
            frame = double(data.meSubset(:, :, i));
            roiMotion(i) = mean(frame(roiMask));
        end

        data.roiMotionEnergySubset{roiIdx} = roiMotion;

        if roiIdx == 1
            cla(ax3);
            hold(ax3, 'on');
        end

        timeVector = (1:length(roiMotion)) / data.frameRate;
        plot(ax3, timeVector, roiMotion, ...
            'Color', data.roiColors{roiIdx}, 'LineWidth', 1.5, ...
            'DisplayName', data.roiNames{roiIdx});

        if roiIdx == 1
            currentTime = data.currentFrame / data.frameRate;
            data.roiTimeLines = xline(ax3, currentTime, 'r', 'LineWidth', 2);
            windowSize = min(60, timeVector(end));
            xlim(ax3, [0, windowSize]);
        end

        legend(ax3, 'Location', 'best');
    end

    function displayFrame(frameNum, updateROIZoom)
        if isempty(data.videoSubset)
            return;
        end

        rawFrame = data.videoSubset(:, :, frameNum);
        if ~isempty(data.rawImageHandle) && isvalid(data.rawImageHandle)
            data.rawImageHandle.CData = rawFrame;
        else
            data.rawImageHandle = imagesc(axRaw, rawFrame);
            colorbar(axRaw);
            axis(axRaw, 'image');
            colormap(axRaw, 'gray');
        end
        title(axRaw, sprintf('Raw Video - Frame %d', frameNum));

        if frameNum <= size(data.meSubset, 3)
            meFrame = data.meSubset(:, :, frameNum);
            if ~isempty(data.meImageHandle) && isvalid(data.meImageHandle)
                data.meImageHandle.CData = meFrame;
            else
                data.meImageHandle = imagesc(axME, meFrame);
                colorbar(axME);
                axis(axME, 'image');
                colormap(axME, 'hot');
            end
            title(axME, sprintf('Motion Energy - Frame %d', frameNum));
        end

        currentTime = frameNum / data.frameRate;
        if ~isempty(data.timeLine) && isvalid(data.timeLine)
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
        timeVector = (1:size(data.meSubset,3)) / data.frameRate;
        if currentTime > windowSize/2 && currentTime < timeVector(end) - windowSize/2
            xlim(ax2, [currentTime - windowSize/2, currentTime + windowSize/2]);
            xlim(ax3, [currentTime - windowSize/2, currentTime + windowSize/2]);
        end

        frameLabel.Text = sprintf('Frame: %d/%d (%.2fs)', ...
            frameNum, data.subsetFrames, currentTime);

        if updateROIZoom
            drawROIsOnAxes();
            updateZoomedROI(frameNum);
        end
    end

    function drawROIsOnAxes()
        warnState = warning('off', 'all');

        hold(axRaw, 'on');
        hold(axProj, 'on');

        for i = 1:length(data.rois)
            if isvalid(data.rois{i})
                try
                    data.rois{i}.Color = data.roiColors{i};
                    data.rois{i}.LineWidth = 2;

                    rect = data.rois{i}.Position;
                    rectangle(axProj, 'Position', rect, ...
                        'EdgeColor', data.roiColors{i}, 'LineWidth', 2);
                catch
                end
            end
        end

        hold(axRaw, 'off');
        hold(axProj, 'off');

        warning(warnState);
    end

    function updateZoomedROI(frameNum)
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

        frame = data.videoSubset(:, :, frameNum);

        rect = roi.Position;
        xStart = max(1, floor(rect(1)));
        xEnd = min(size(frame,2), ceil(rect(1) + rect(3)));
        yStart = max(1, floor(rect(2)));
        yEnd = min(size(frame,1), ceil(rect(2) + rect(4)));

        zoomedFrame = frame(yStart:yEnd, xStart:xEnd);

        if ~isempty(data.zoomImageHandle) && isvalid(data.zoomImageHandle)
            data.zoomImageHandle.CData = zoomedFrame;
        else
            data.zoomImageHandle = imagesc(axZoom, zoomedFrame);
            colorbar(axZoom);
            axis(axZoom, 'image');
            colormap(axZoom, 'gray');
        end

        title(axZoom, sprintf('%s Zoom - Frame %d', data.roiNames{data.selectedROI}, frameNum));

        hold(axZoom, 'on');
        rectangle(axZoom, 'Position', [1, 1, size(zoomedFrame,2), size(zoomedFrame,1)], ...
            'EdgeColor', data.roiColors{data.selectedROI}, 'LineWidth', 2);
        hold(axZoom, 'off');
    end

    function frameSliderCallback(src, ~)
        data.currentFrame = round(src.Value);
        displayFrame(data.currentFrame, true);
        drawnow;
    end

    function playCallback(~, ~)
        data.isPlaying = true;
        playBtn.Enable = 'off';
        stopBtn.Enable = 'on';

        targetFPS = data.frameRate;
        frameInterval = 1 / targetFPS;
        displayEveryNFrames = max(1, round(targetFPS / 30));

        frameCount = 0;
        lastTime = tic;

        while data.isPlaying && data.currentFrame <= data.subsetFrames
            if toc(lastTime) >= frameInterval
                if mod(frameCount, displayEveryNFrames) == 0
                    displayFrame(data.currentFrame, true);
                    frameSlider.Value = data.currentFrame;
                end

                data.currentFrame = data.currentFrame + 1;
                frameCount = frameCount + 1;
                lastTime = tic;
            end

            drawnow limitrate;

            if ~data.isPlaying
                break;
            end
        end

        data.isPlaying = false;
        playBtn.Enable = 'on';
        stopBtn.Enable = 'on';

        if data.currentFrame > data.subsetFrames
            data.currentFrame = 1;
        end

        displayFrame(data.currentFrame, true);
        drawnow;
    end

    function stopCallback(~, ~)
        data.isPlaying = false;
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
        data.roiRects(idx) = [];
        data.roiMotionEnergySubset(idx) = [];

        roiListBox.Items = data.roiNames;
        updateROITable();

        data.selectedROI = [];
        data.selectedROIIdx = [];

        cla(ax3);
        hold(ax3, 'on');
        for i = 1:length(data.roiMotionEnergySubset)
            if ~isempty(data.roiMotionEnergySubset{i})
                timeVector = (1:length(data.roiMotionEnergySubset{i})) / data.frameRate;
                plot(ax3, timeVector, data.roiMotionEnergySubset{i}, ...
                    'Color', data.roiColors{i}, 'LineWidth', 1.5, ...
                    'DisplayName', data.roiNames{i});
            end
        end
        legend(ax3, 'Location', 'best');
        hold(ax3, 'off');

        displayFrame(data.currentFrame, true);
    end

    function clearAllROIsCallback(~, ~)
        for i = 1:length(data.rois)
            if isvalid(data.rois{i})
                delete(data.rois{i});
            end
        end

        data.rois = {};
        data.roiNames = {};
        data.roiColors = {};
        data.roiRects = {};
        data.roiMotionEnergySubset = {};
        data.selectedROI = [];
        data.selectedROIIdx = [];

        roiListBox.Items = {};
        updateROITable();

        cla(axZoom);
        cla(ax3);
        data.zoomImageHandle = [];
        data.roiTimeLines = [];

        deleteROIBtn.Enable = 'off';
        clearROIBtn.Enable = 'off';

        displayFrame(data.currentFrame, true);
    end

    function roiSelectionCallback(src, ~)
        selectedName = src.Value;
        if isempty(selectedName)
            return;
        end

        data.selectedROI = find(strcmp(data.roiNames, selectedName), 1);
        data.selectedROIIdx = data.selectedROI;
        displayFrame(data.currentFrame, true);
    end

    function roiTableSelectionCallback(~, event)
        if ~isempty(event.Indices)
            data.selectedROIIdx = event.Indices(1);
            data.selectedROI = data.selectedROIIdx;

            if data.selectedROIIdx <= length(data.roiNames)
                roiListBox.Value = data.roiNames{data.selectedROIIdx};
                displayFrame(data.currentFrame, true);
            end
        end
    end

    function updateROITable()
        if isempty(data.rois)
            roiTable.Data = {};
            return;
        end

        tableData = cell(length(data.rois), 2);
        for i = 1:length(data.rois)
            tableData{i, 1} = data.roiNames{i};
            tableData{i, 2} = sprintf('[%.1f %.1f %.1f]', data.roiColors{i});
        end
        roiTable.Data = tableData;
    end

    function processVideoCallback(~, ~)
        if isempty(data.outputDir)
            uialert(fig, 'Please set output directory first.', 'No Output Dir');
            return;
        end

        processFullVideoSVD();
    end

    function processFullVideoSVD()
        [~, videoName, ~] = fileparts(data.videoName);

        try
            progressLabel.Text = sprintf('Processing %s with SVD...', videoName);
            drawnow;

            outputPath = fullfile(data.outputDir, [videoName '_motionSVD.mat']);

            vidObj = VideoReader(data.videoPath);
            totalFrames = vidObj.NumFrames;
            frameRate = vidObj.FrameRate;
            frameHeight = vidObj.Height;
            frameWidth = vidObj.Width;
            nPixels = frameHeight * frameWidth;

            %% Step 1: Compute average motion energy
            progressLabel.Text = 'Computing average motion energy...';
            progressPercent.Text = '5%';
            drawnow;

            nFramesAvg = min(1000, totalFrames);
            avgMotSum = zeros(frameHeight, frameWidth);

            vidObj.CurrentTime = 0;
            prevFrame = [];

            for i = 1:nFramesAvg
                if hasFrame(vidObj)
                    frame = readFrame(vidObj);
                    if size(frame, 3) == 3
                        frame = rgb2gray(frame);
                    end

                    if ~isempty(prevFrame)
                        motFrame = abs(int16(frame) - int16(prevFrame));
                        avgMotSum = avgMotSum + double(motFrame);
                    end
                    prevFrame = frame;
                end
            end

            avgmot = avgMotSum / (nFramesAvg - 1);
            avgmot_vec = avgmot(:);

            %% Step 2: Compute whole-frame SVD motion masks
            progressLabel.Text = 'Computing whole-frame SVD...';
            progressPercent.Text = '10%';
            drawnow;

            chunkSize = 500;
            nChunksForSVD = min(10, floor(totalFrames / chunkSize));

            uMot = [];

            for j = 1:nChunksForSVD
                startFrame = (j-1)*chunkSize + 1;
                endFrame = min(j*chunkSize + 1, totalFrames);
                nFramesChunk = endFrame - startFrame + 1;

                progressPercent.Text = sprintf('%.0f%%', 10 + (j/nChunksForSVD)*15);
                drawnow;

                vidObj.CurrentTime = (startFrame-1) / frameRate;
                F = zeros(frameHeight, frameWidth, nFramesChunk, 'uint8');

                for i = 1:nFramesChunk
                    if hasFrame(vidObj)
                        frame = readFrame(vidObj);
                        if size(frame, 3) == 3
                            frame = rgb2gray(frame);
                        end
                        F(:, :, i) = frame;
                    end
                end

                M = zeros(nPixels, nFramesChunk-1);
                for i = 2:nFramesChunk
                    motFrame = abs(int16(F(:,:,i)) - int16(F(:,:,i-1)));
                    M(:, i-1) = double(motFrame(:));
                end

                M = M - avgmot_vec;

                [u,~,~] = svd(M, 'econ');
                uMot = [uMot, u];

                clear F M;
            end

            [uMot,~,~] = svd(uMot, 'econ');
            uMotMask = normc(uMot(:, 1:min(500, size(uMot,2))));

            clear uMot;

            %% Step 3: Project full video onto motion masks
            progressLabel.Text = 'Projecting video onto motion masks...';
            drawnow;

            nChunks = ceil((totalFrames-1) / chunkSize);
            motSVD = zeros(totalFrames-1, size(uMotMask, 2));

            vidObj.CurrentTime = 0;
            frameIdx = 1;

            for j = 1:nChunks
                startFrame = (j-1)*chunkSize + 1;
                endFrame = min(j*chunkSize + 1, totalFrames);
                nFramesChunk = endFrame - startFrame + 1;

                progressLabel.Text = sprintf('Projecting chunk %d/%d...', j, nChunks);
                progressPercent.Text = sprintf('%.0f%%', 25 + (j/nChunks)*35);
                drawnow;

                vidObj.CurrentTime = (startFrame-1) / frameRate;
                F = zeros(frameHeight, frameWidth, nFramesChunk, 'uint8');

                for i = 1:nFramesChunk
                    if hasFrame(vidObj)
                        frame = readFrame(vidObj);
                        if size(frame, 3) == 3
                            frame = rgb2gray(frame);
                        end
                        F(:, :, i) = frame;
                    end
                end

                M = zeros(nPixels, nFramesChunk-1);
                for i = 2:nFramesChunk
                    motFrame = abs(int16(F(:,:,i)) - int16(F(:,:,i-1)));
                    M(:, i-1) = double(motFrame(:));
                end

                motSVD0 = (M - avgmot_vec)' * uMotMask;

                endIdx = frameIdx + size(motSVD0, 1) - 1;
                motSVD(frameIdx:endIdx, :) = motSVD0;
                frameIdx = endIdx + 1;

                clear F M motSVD0;
            end

            %% Step 4: Load full motion data for ROI SVD
            progressLabel.Text = 'Loading motion data for ROI analysis...';
            progressPercent.Text = '60%';
            drawnow;

            motionData = zeros(frameHeight, frameWidth, totalFrames-1, 'single');

            vidObj.CurrentTime = 0;
            prevFrame = [];
            idx = 1;

            for i = 1:totalFrames
                if hasFrame(vidObj)
                    frame = readFrame(vidObj);
                    if size(frame, 3) == 3
                        frame = rgb2gray(frame);
                    end

                    if ~isempty(prevFrame)
                        motionData(:,:,idx) = single(abs(int16(frame) - int16(prevFrame)));
                        idx = idx + 1;
                    end
                    prevFrame = frame;
                end

                if mod(i, 500) == 0
                    progressPercent.Text = sprintf('%.0f%%', 60 + (i/totalFrames)*15);
                    drawnow limitrate;
                end
            end

            motionData = motionData - single(avgmot);

            %% Step 5: Per-ROI SVD and traces
            progressLabel.Text = 'Computing per-ROI SVD...';
            progressPercent.Text = '75%';
            drawnow;

            bodyPartData = struct();

            if ~isempty(data.rois)
                motionFlat = reshape(motionData, [], size(motionData,3));

                for roiIdx = 1:length(data.rois)
                    rect = data.roiRects{roiIdx};
                    roiMask = false(frameHeight, frameWidth);
                    xStart = max(1, round(rect(1)));
                    xEnd = min(frameWidth, round(rect(1) + rect(3)));
                    yStart = max(1, round(rect(2)));
                    yEnd = min(frameHeight, round(rect(2) + rect(4)));
                    roiMask(yStart:yEnd, xStart:xEnd) = true;

                    idx_roi = find(roiMask(:));
                    F_ROI = motionFlat(idx_roi, :);
                    F_ROI = F_ROI - mean(F_ROI, 2);

                    [~,S,~] = svd(F_ROI, 'econ');
                    singvals = diag(S);

                    roiMotion = mean(F_ROI, 1)';

                    roiName = data.roiNames{roiIdx};
                    bodyPartData.([roiName '_Position']) = rect;
                    bodyPartData.([roiName '_mask']) = roiMask;
                    bodyPartData.([roiName '_motion']) = roiMotion;
                    bodyPartData.([roiName '_singvals']) = singvals;

                    progressPercent.Text = sprintf('%.0f%%', 75 + (roiIdx/length(data.rois))*20);
                    drawnow;
                end
            end

            clear motionData motionFlat;

            %% Save output
            progressLabel.Text = 'Saving results...';
            progressPercent.Text = '95%';
            drawnow;

            saveData = struct();
            saveData.motSVD = motSVD;
            saveData.uMotMask = uMotMask;
            saveData.avgmot = avgmot;
            saveData.frameRate = frameRate;
            saveData.totalFrames = totalFrames;
            saveData.timeVector = (1:totalFrames-1)' / frameRate;
            saveData.frameSize = [frameHeight, frameWidth];
            saveData.bodyPartData = bodyPartData;

            save(outputPath, '-struct', 'saveData', '-v7.3');

            progressLabel.Text = 'Processing complete!';
            progressPercent.Text = '100%';
            statusLabel.Text = sprintf('Saved: %s', videoName);

            uialert(fig, sprintf('Processing complete!\n\nSaved to: %s\n\nContents:\n- motSVD: [%d x %d] whole-frame\n- uMotMask: [%d x %d]\n- Body parts: %d ROIs with SVD\n- Per-ROI singular values saved', ...
                outputPath, size(motSVD,1), size(motSVD,2), ...
                size(uMotMask,1), size(uMotMask,2), length(data.rois)), ...
                'Success', 'Icon', 'success');

        catch ME
            uialert(fig, sprintf('Error processing video: %s', ME.message), 'Processing Error');
            progressLabel.Text = 'Processing failed';
            progressPercent.Text = '';
        end
    end

end
