function [idx1,idx2,idx3] = powerSpecStat(pxSpecs)
thetaP = pxSpecs.broadband((pxSpecs.f>4 & pxSpecs.f<12),:);
betaP = pxSpecs.broadband((pxSpecs.f>12 & pxSpecs.f<32),:);
gammaP1 = pxSpecs.broadband((pxSpecs.f>32 & pxSpecs.f<58),:); % Cutout 60 Hz
gammaP2 =  pxSpecs.broadband((pxSpecs.f>62 & pxSpecs.f<100),:);
gammaP = [gammaP1;gammaP2];
idx1 = mean(thetaP,1);idx2 = mean(betaP,1);idx3 = mean(gammaP,1);
total = [mean(idx1) mean(idx2) mean(idx3)];
err = [std(idx1)/sqrt(length(idx1)) std(idx2)/sqrt(length(idx2)) std(idx3)/sqrt(length(idx3))];
figure,bar(total),hold on
errorbar(1:3,total,err)
idx1 = idx1';idx2 = idx2';idx3 = idx3';