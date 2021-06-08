function [sortedPlaceFields,normPlaceFields] = placefieldSort(allplaceFields)

placeFields = allplaceFields;
FRmax = max(placeFields,[],2);
count = 1;
for i = 1:size(placeFields,1)
    if isempty(find(placeFields(i,:)==FRmax(i,1)))
        idx(i,1) = i; % Index of neuron
        idx(i,2) = 0; % Position in space
        count = count+1;
    elseif FRmax(i,1) == 0
        count = count+1;
    else
        idx(count,1) = i;
        idx(count,2) = find(placeFields(i,:)==FRmax(i,1));
        count = count+1;
    end
end

% Now to sort the idx from least to greatest
[~,idx_sort] = sort(idx(:,2),1); % Ascending list of locations
sortedPlaceFields = placeFields(idx_sort,:);
padding = length(find(FRmax==0));
sortedPlaceFields = sortedPlaceFields(padding+1:end,:);

% Normalize Place Fields for visualization
for ii = 1:size(sortedPlaceFields,1)
    maxV = max(sortedPlaceFields(ii,:));
    minV = min(sortedPlaceFields(ii,:));
    normPlaceFields(ii,:) = (sortedPlaceFields(ii,:)-minV)./(maxV-minV);
    % Clean up noise
end



   