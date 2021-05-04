function linescan = spineScan(fn)
if nargin < 1 || strcmp(fn,'')
    [filename,pathname] = uigetfile('*.csv');
    fn = [pathname filename];
else 
    filename=fn;
end

T = table2array(readtable(fn));
T(:,1) = [];
for i = 1:size(T,2)
    pow(:,i) = fastrms(T(:,i),5);
    mLine(i) = mean(pow(:,i));
%     meLine(i) = median(pow(i));
%     mmLineRatio(i) = mLine(i)./meLine(i);
end
% mmLineRatio(mLine<1) = 0;
% mmLineRatio(meLine<1) = 0;

[maxVal,loc] = max(mLine);
disp([ 'Best linescan detected: ' num2str(loc)]);
linescan = T(:,loc);
figure()
plot(linescan,'LineWidth',1); box off, axis tight;
end



