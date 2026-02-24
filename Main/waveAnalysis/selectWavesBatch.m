function [WavesOut] = selectWavesBatch(Waves,stEvalP,spEvalP)

WavesOut = Waves;

for i=1:size(Waves,2)
    if isempty(Waves(i).evaluationPoints)
        WavesOut(i).wavePresent = WavesOut(i).wavePresent(stEvalP:spEvalP);
        WavesOut(i).evaluationPoints = [];
        WavesOut(i).PGD = WavesOut(i).PGD(stEvalP:spEvalP);
        WavesOut(i).rho = [];
        WavesOut(i).source = [];
        WavesOut(i).speed = [];
        WavesOut(i).waveDir = [];
        WavesOut(i).wavelength = [];
    else
        posDel = find(Waves(i).evaluationPoints < stEvalP | Waves(i).evaluationPoints > spEvalP);
        WavesOut(i).wavePresent = WavesOut(i).wavePresent(stEvalP:spEvalP);
        WavesOut(i).evaluationPoints(posDel) = [];
        WavesOut(i).PGD = WavesOut(i).PGD(stEvalP:spEvalP);
        WavesOut(i).rho(posDel) = [];
        WavesOut(i).source(posDel,:) = [];
        WavesOut(i).speed(posDel) = [];
        WavesOut(i).waveDir(posDel) = [];
        WavesOut(i).wavelength(posDel) = [];
    end
end

