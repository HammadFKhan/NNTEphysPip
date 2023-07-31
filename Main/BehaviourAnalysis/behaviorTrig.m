function VmTrigAlign = behaviorTrig(out,VmTotal)
%% Behavior Trig
% Fix DC offset
Trig = out.S{5};
TrigTotal = [];
for i = 1:size(Trig,2)
    TrigTotal = [TrigTotal;Trig(:,i)];
end

flag = find(TrigTotal>2); % Voltage swing from ~0 to ~4 mV;
startflag = flag(1); % grab first indice
endflag = flag(end); %grab last indice
VmTrigAlign = VmTotal(startflag:endflag); % Only include trace when trig starts


%%
Frame = out.S{3};
FrameTotal = [];
for i = 1:size(Frame,2)
    FrameTotal = [FrameTotal;Frame(:,i)];
end


end

