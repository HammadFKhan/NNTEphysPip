function Waves = zscoreWavesSpeedPGD(Waves, parameters)
% Z-score speed and PGD for Waves.* sub-structs.
% Baseline mean/std are computed from waves between 0 and 1500 ms
% using selectWaves, then stored as new fields .zSpeed and .zPGD.

    % Time window in samples for baseline (0–1500 ms relative to cue)
    srt = 1;
    stp = round(parameters.windowBeforeCue * parameters.Fs);

    % Helper: z-score one field across an array of wave-trial structs
    function wavesOut = zscoreField(wavesIn, fieldName, newFieldName)
        wavesOut = wavesIn;
        if isempty(wavesIn)
            return;
        end

        % 1) Collect baseline values from 0–1500 ms using selectWaves
        wavesBaseline = selectWaves(wavesIn, srt, stp);

        baselineVals = [];
        for kk = 1:numel(wavesBaseline)
            if isfield(wavesBaseline(kk), fieldName) && ~isempty(wavesBaseline(kk).(fieldName))
                vals = wavesBaseline(kk).(fieldName);
                baselineVals = [baselineVals, vals(:).'];
            end
        end

        if isempty(baselineVals)
            warning('No baseline values found for field "%s"; skipping z-scoring.', fieldName);
            return;
        end

        mu    = mean(baselineVals, 'omitnan');
        sigma = std(baselineVals, 0, 'omitnan');
        if sigma == 0 || isnan(sigma)
            warning('Zero/NaN std for field "%s"; skipping z-scoring.', fieldName);
            return;
        end

        % 2) Apply z-score to every wave in the full array
        for kk = 1:numel(wavesOut)
            if isfield(wavesOut(kk), fieldName) && ~isempty(wavesOut(kk).(fieldName))
                wavesOut(kk).(newFieldName) = (wavesOut(kk).(fieldName) - mu) ./ sigma;
            else
                wavesOut(kk).(newFieldName) = [];
            end
        end
    end

    % Process each condition field if present
    condFields = {'wavesHit','wavesMiss','wavesFA','wavesHitReward','wavesMIHit','wavesMIFA','wavesOptoCueHit','wavesOptoCueMiss'};

    for c = 1:numel(condFields)
        fld = condFields{c};
        if isfield(Waves, fld) && ~isempty(Waves.(fld))
            Waves.(fld) = ensureWaveStart(Waves.(fld));
            Waves.(fld) = zscoreField(Waves.(fld), 'speed', 'zSpeed');
            Waves.(fld) = zscoreField(Waves.(fld), 'PGD',   'zPGD');
        end
    end
end
