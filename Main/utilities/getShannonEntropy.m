function [H] = getShannonEntropy(x,nBins,minVal,maxVal)

binEdges = linspace(minVal,maxVal,nBins+1);

% Binning the array 
binnedX = discretize(x,binEdges);

pk = zeros(nBins,1);

Hk = zeros(nBins,1);

for k=1:nBins
    pk(k) = (sum(binnedX==k))/size(x,2);
    if pk(k) == 0 
        Hk(k) = 0;
    else
        Hk(k) = -1*pk(k)*(log2(pk(k)));
    end
end

H = sum(Hk);

end
