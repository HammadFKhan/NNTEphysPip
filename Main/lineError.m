function lineError(x,y,varargin)
if strcmp(varargin,'std')
        errorType = 0;
elseif strcmp(varargin,'ste')
        errorType = 1;
else
    warning('No error type given, defaulting to std')
    errorType = 1;     
end

if errorType == 0
    error = std(y);
    yM = mean(y,1);
%     x = x';
    curve1 = yM + error;
    curve2 = yM - error;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, [0.75 0.75 0.75]);
    hold on;
    plot(x, yM, 'k', 'LineWidth', 2);
else
    yM = mean(y,1);
%     x = x';
    error = std(y);
    error = error/sqrt(length(yM));
    curve1 = yM + error;
    curve2 = yM - error;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, [0.75 0.75 0.75]);
    hold on;
    plot(x, yM, 'k', 'LineWidth', 2);
end
end