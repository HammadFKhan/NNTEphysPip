function plotPlacefield(Spikes)

%Parse
placeFieldMap = Spikes.PlaceFields.placeFieldMap;
placeFieldperTrial = Spikes.PlaceFields.placeFieldperTrial;
allplaceFields = Spikes.PlaceFields.allplaceFields;
sortedPlaceFields = Spikes.PlaceFields.sortedPlaceFields;
normPlaceFields = Spikes.PlaceFields.normPlaceFields;
populationVector = Spikes.PlaceFields.populationVector;
normplaceFieldperTrial = Spikes.PlaceFields.normplaceFieldperTrial;
avg_centeredPlaceFields = Spikes.PlaceFields.avg_centeredPlaceFields;
trackLength = 200;
disp(['Track Length set to ' num2str(trackLength) ' cm'])
pattern = 70;
disp(['New Pattern set at ' num2str(pattern) ' cm'])
figure('Name','Place Field Candidates');

for neuron = 1:size(Spikes.PlaceFields.placeField,2)
    position = 1:trackLength;
    s = ceil(sqrt(size(Spikes.PlaceFields.placeField,2)));
    subplot(s,s,neuron),spike_map(placeFieldMap(:,:,neuron),position);
    axis tight, axis off, colorbar off;
    title(['PF ' num2str(neuron)]);
end
% Place Field Map per Trial
[grad,~]=colorGradient([1 1 1],[.1 .1 .1],64); %color map gradient for nice visuals
figure('Name','Place Field per Trial')
count = 1;
for perTrial = 17:20
    subplot(1,4,count),spike_map(normplaceFieldperTrial{perTrial}(:,1:100),1:100,grad);
    hold on,colorbar( 'off' ),set(gca,'ytick',[]),...
set(gca,'yticklabel',[]),xlabel(''),ylabel('')
    title(['Trial ' num2str(count)])
    count = count+1;
end

figure('name','PF Curves')
for trial = [7:8]
    for neuron = 1:size(normplaceFieldperTrial{1,trial},1)
    plot(normplaceFieldperTrial{1,trial}(neuron,1:100),'LineWidth',2), hold on,box off;
    end
end
ylim([0 4]);
% Plot Field Width
figure('Name','Centered Place Field Width')
count = 1;
% for perTrial = [9,12,10]
%     plot(Smooth(avg_centeredPlaceFields(perTrial,1:40),5),'LineWidth',3);
%     hold on,box off,colorbar( 'off' ),set(gca,'ytick',[]),...
% set(gca,'yticklabel',[]),set(gca,'xtick',[]),...
% set(gca,'xticklabel',[])
% end

% figure,plot(Smooth(mean(avg_centeredPlaceFields([14 15 16],1:40),1),10))
% figure, subplot(4,1,[1 3]),spike_map(placeFieldMap(9:16,:,6),position);


figure('Name','All Place Fields across Trials')
spike_map(allplaceFields,1:size(allplaceFields,2));

% Sort Place Fields
figure('Name','Sorted Place Fields'),spike_map(sortedPlaceFields,position);
figure('Name','Normalized Place Fields'),spike_map(normPlaceFields(:,1:100),1:100);hold on;
plot(ones(size(normPlaceFields,1),1)*pattern, 1:size(normPlaceFields,1),'--w','LineWidth',2);
plot(ones(size(normPlaceFields,1),1)*pattern+20, 1:size(normPlaceFields,1),'--w','LineWidth',2); hold off
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/PlaceFieldmap.eps', '-r250');
% Population Vector
figure('Name','Spatial Scale Factor'),htmp(populationVector,20);

end
