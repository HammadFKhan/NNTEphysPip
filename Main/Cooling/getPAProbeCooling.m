function [BaselinePA,CooledPA] = getPAProbeCooling(LFP,IntanBehaviour,z_score,nPerm,plotFlag,parameters)

% Reference : Spontaneous travelling cortical waves gate perception in 
% behaving primates, Nature, 2020 

% Similair to getPA but adjusted for silicon probe data. (HK 06/2024)
nElectrodes = parameters.rows*parameters.cols;

% Make cooled structure and then we just call getPAProbe twice for each
% temperature type. Easy peasy

LFP.Cooled.hitxgp = LFP.hitxgp(IntanBehaviour.hitTemp<-10);
LFP.Cooled.missxgp = LFP.missxgp(IntanBehaviour.missTemp<-10);
LFP.Cooled.MIhitxgp = LFP.MIhitxgp(IntanBehaviour.hitTemp<-10);
LFP.Cooled.MIFAxgp = LFP.MIFAxgp(IntanBehaviour.FATemp<-10);

LFP.Baseline.hitxgp = LFP.hitxgp(IntanBehaviour.hitTemp>-10);
LFP.Baseline.missxgp = LFP.missxgp(IntanBehaviour.missTemp>-10);
LFP.Baseline.MIhitxgp = LFP.MIhitxgp(IntanBehaviour.hitTemp>-10);
LFP.Baseline.MIFAxgp = LFP.MIFAxgp(IntanBehaviour.FATemp>-10);

% Then we just call the function twice here, jjeje
[BaselinePA] = getPAProbe(LFP.Baseline,IntanBehaviour,z_score,nPerm,plotFlag,parameters);
[CooledPA] = getPAProbe(LFP.Cooled,IntanBehaviour,z_score,nPerm,plotFlag,parameters);

end