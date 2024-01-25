function [PPL] = calPPL(xgp,parameters)

% Reference : Propagating waves mediate information transfer in the motor
% cortex, Nature Neuro, 2006 
   
totaltime = parameters.windowAfterCue + parameters.windowBeforeCue + parameters.ts;
time = [parameters.ts:parameters.ts:totaltime]; 

nTrials = size(xgp,2);
% allXGP = zeros(nTrials ,size(time,2));
% for trialno=1:nTrials 
%     allXGP(trialno,:) = angle(xgp{1,trialno});
% end

% N = ceil(exp(0.626 + (0.4*log(nTrials-1)))); % number of bins
N = 18;
% Ref - Comparison of Hilbert transform and wavelet methods for the
%       analysis of neuronal synchrony, JNeuroMethods, 2001

Hmax = log2(N);

H = zeros(size(xgp{1,1},1),size(xgp{1,1},2),size(xgp{1,1},3));
PPL = zeros(size(xgp{1,1},1),size(xgp{1,1},2),size(xgp{1,1},3));

for t=1:size(xgp{1,1},3)
    for i=1:size(xgp{1,1},1)
        for j=1:size(xgp{1,1},2)
            x = cellfun(@(s) angle(s(i,j,t)),xgp);
            if sum(isnan(x))>0
                PLL(i,j,t) = NaN;
                H(i,j,t) = NaN;
            else
                H(i,j,t) = getShannonEntropy(x,N,-pi,pi);
                PPL(i,j,t) = 100*(1-(H(i,j,t)/Hmax));
            end
        end
    end
end

