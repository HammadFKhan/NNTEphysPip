function Waves = ensureWaveStart(Waves)

for ii = 1:numel(Waves)
    % If waveStart already exists, skip
    if isfield(Waves, 'waveStart') && ~isempty(Waves(ii).waveStart)
        continue;
    end

    % If no wavePresent field, nothing to do
    if ~isfield(Waves, 'wavePresent') || isempty(Waves(ii).wavePresent)
        continue;
    end

    % Create waveStart from wavePresent
    wp = Waves(ii).wavePresent(:).';          % row vector
    ws = zeros(size(wp));                     % init

    if any(wp)
        % wave start = first element that is 1, or 1 preceded by 0
        ws(wp == 1 & [true, wp(1:end-1) == 0]) = 1;
    end

    Waves(ii).waveStart = ws;
end
end
