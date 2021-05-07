function T = spineScan(fn)
if nargin < 1 || strcmp(fn,'')
    [filename,pathname] = uigetfile('*.csv');
    fn = [pathname filename];
else 
    filename=fn;
end

T = table2array(readtable(fn));
T(:,1) = [];
prev = T(1:3000,:);
meanLine = prev-mean(prev);

figure('Name','Line Channels')
num = ceil(size(T,2)/3);
for i = 1:size(T,2)
    normLine(:,i) = (meanLine(:,i)-min(meanLine,[],'all'))/((max(meanLine,[],'all'))-min(meanLine,[],'all'));
%     pow(:,i) = rms(T(:,i));
%     mLine(i) = mean(pow(:,i));
%     meLine(i) = median(pow(i));
%     mmLineRatio(i) = mLine(i)./meLine(i);
%     
%     powC(:,i) = rms(normLine(:,i));
%     mLineC(i) = mean(powC(:,i));
%     meLineC(i) = median(powC(i));
%     mmLineRatioC(i) = mLineC(i)./meLineC(i);
    subplot(3,num,i),plot(smooth(normLine(:,i),0.01,'sgolay'),'LineWidth',1); 
    box off, axis tight; title(['Line Channel ' num2str(i)]);
    set(gca,'YTick','','YTickLabel','');
end
% mmLineRatio(mLine<1) = 0;
% mmLineRatio(meLine<1) = 0;

% [maxVal,loc] = max(mLineC);
% disp([ 'Best linescan detected: ' num2str(loc)]);
% linescan = T(:,10);
% figure()
% plot(linescan,'LineWidth',1); box off, axis tight;
end



