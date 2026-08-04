function [neuralTrajSim,neuralTrajdiff,rprimehnorm,rprimemnorm] = neuralTrajDiff(r1,r2,varargin)
if strcmp(varargin,'initial'),initCalc = 1;else, initCalc = 0;end
if initCalc %r'(t)/||r'(t)|| * r(0)-ri(t)/||r(0)-ri(t)||
    rprimeh = diff(r1);
    rprimehnorm = rprimeh/norm(rprimeh);
    rprimem = r2(1,:)-r1; % control condition for reference
    rprimemnorm = rprimem/norm(rprimem);
    neuralTrajSim = dot(rprimehnorm',rprimemnorm(1:end-1,:)');
    neuralTrajdiff = rprimeh-rprimem;
else
    rprimeh = diff(r1);
    rprimehnorm = rprimeh/norm(rprimeh);
    rprimem = diff(r2);
    rprimemnorm = rprimem/norm(rprimem);
    neuralTrajSim = dot(rprimehnorm',rprimemnorm');
    neuralTrajdiff = rprimeh-rprimem;
end
end
