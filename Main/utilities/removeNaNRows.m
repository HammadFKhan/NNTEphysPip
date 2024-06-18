function newArray = removeNaNRows(array)
    % Find rows containing NaN values
    nanRows = any(isnan(array), 2);
    
    % Remove rows with NaN values
    newArray = array(~nanRows, :);
end
