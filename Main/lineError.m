function lineError(x,y,varargin)
if strcmp(varargin,'std')
        errorType = 0;
elseif strcmp(varargin,'ste')
        errorType = 1;
else
    warning('No error type given, defaulting to std')
    errorType = 0;     
end

if errorType == 0
    y = mean(y,1);
    error = std(y);
    curve1 = y + error;
    curve2 = y - error;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, [0.75 0.75 0.75]);
    hold on;
    plot(x, y, 'k', 'LineWidth', 2);
else
    y = mean(y,1);
    error = std(y);
    error = error/sqrt(length(y));
    curve1 = y + error;
    curve2 = y - error;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, [0.75 0.75 0.75]);
    hold on;
    plot(x, y, 'k', 'LineWidth', 2);
end
end