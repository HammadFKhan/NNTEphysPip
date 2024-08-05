% AP scaling based on Ra
function dataScale = APscale(data)
idx = find(data>0);
dataScale = data;
scaleFactor = 60/max(data(idx));
dataScale(idx) = dataScale(idx)*scaleFactor;
figure,
subplot(2,1,1),plot(data),box off,set(gca,'Tickdir','out'),axis tight
subplot(2,1,2),plot(dataScale,'r'),box off,set(gca,'Tickdir','out'),axis tight,linkaxes