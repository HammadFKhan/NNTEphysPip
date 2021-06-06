function sortedPlaceFields = placefieldSort(allplaceFields)

placeFields = allplaceFields;
FRmax = max(placeFields,[],2);
for i = 1:size(placeFields,1)
    if isempty(find(placeFields(i,:)==FRmax(i,1)))
        idx(i,1) = i; % Index of neuron
        idx(i,2) = 0; % Position in space
    else
        idx(i,1) = i;
        idx(i,2) = find(placeFields(i,:)==FRmax(i,1));
    end
end

% Now to sort the idx from least to greatest
[~,idx_sort] = sort(idx(:,2),1); % Ascending list of locations
sortedPlaceFields = placeFields(idx_sort,:);
sortedPlaceFields = sortedPlaceFields;


   