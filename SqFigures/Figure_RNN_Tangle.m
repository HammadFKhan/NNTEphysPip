load('C:\Users\khan332\Documents\GitHub\ssm\sequenceProject\RNN_model_data.mat')
load('C:\Users\khan332\Documents\GitHub\ssm\sequenceProject\RNN_model_data_lowTangle.mat')
 %% Compute tangle for simulated trajectories
 % Tangle input as X: latent_dim x time x trials
X = permute(RNN_output_tangled.target,[3,2,1]);
Qs = computeTangling(X,0.0001);
Qs = max(Qs,[],1);
X = permute(RNN_output_tangled.RNN_output,[3,2,1]);
Qrnn = computeTangling(X,0.0001);
Qrnn = max(Qrnn,[],1);
temp = nan(max([length(Qs),length(Qrnn)]),2);
temp(1:length(Qs),1) = Qs;
temp(1:length(Qrnn),2) = Qrnn;
figure,
customBarplot(Qs')
Qs_high = Qs;
 %% Compute tangle for simulated trajectories
 % Tangle input as X: latent_dim x time x trials
X = permute(RNN_output.target,[3,2,1]);
Qs = computeTangling(X,0.0001);
Qs = max(Qs,[],1);
X = permute(RNN_output.RNN_output,[3,2,1]);
Qrnn = computeTangling(X,0.0001);
Qrnn = max(Qrnn,[],1);
temp = nan(max([length(Qs),length(Qrnn)]),2);
temp(1:length(Qs),1) = Qs;
temp(1:length(Qrnn),2) = Qrnn;
figure,
customBarplot([Qs',Qs_high']);ylim([0 4]),set(gca,'tickdir','out'),axis square
%%
X = smoothdata(permute(RNN_output_tangled.target,[3,2,1]),2,'gaussian',20);
figure,plot3(squeeze(X(1,:,:)),squeeze(X(2,:,:)),squeeze(X(3,:,:)),'color',[0 0 0 0.24])
%%
RNNSq= struct();
files = dir(fullfile('C:\Users\khan332\Documents\GitHub\ssm\sequenceProject\ExternalPulse','*.mat'));
for fileNum = 1:length(files)
    disp(['File number: ' num2str(fileNum)])
    load(fullfile(files(fileNum).folder,files(fileNum).name))
    RNNSq(fileNum).name = files(fileNum).name;
    RNNSq(fileNum).target = RNN_data.target;
    RNNSq(fileNum).input = RNN_data.input;
    RNNSq(fileNum).RNN = RNN_data.RNN_output;
    X = permute(RNN_data.target,[3,2,1]);
    X = smoothdata(X,2,'gaussian',20);
    Qs = computeTangling(X,0.0001);
    Qs = max(Qs,[],1);
    RNNSq(fileNum).Qs = Qs;
    Qtotal{fileNum} = Qs;
end
%%

data = vertcat(Qtotal{:})';

% Calculate the mean and standard deviation for each column
mean_data = mean(data, 1);
std_data = std(data, 0, 1)/sqrt(100);

% Define the x-axis values. You'll need to replace these with your actual 'Ramping input slope' values.
x_values = [1:length(mean_data)];
% Create a new figure
figure;
hold on; % This is important to plot multiple things on the same axes

% Plot the shaded error region
% We'll use the 'fill' function for this. We need to define the upper and lower bounds.
upper_bound = mean_data + std_data;
lower_bound = mean_data - std_data;

% Create a polygon for the fill function. We need to close the shape.
x_fill = [x_values, fliplr(x_values)];
y_fill = [upper_bound, fliplr(lower_bound)];

% Plot the filled region with a light gray color
fill(x_fill, y_fill, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plot the main line (mean values) on top of the shaded region
plot(x_values, mean_data, 'k-', 'LineWidth', 2);

% Customize the plot to match the style of the image
% Set axis limits
xlim([0 6]);
% ylim([0.5 0.85]);

% Add labels and a title
xlabel('Input Sequence', 'FontSize', 12);
ylabel('Neural Tangle', 'FontSize', 12);
title('Ramping vs proportion of switching trials', 'FontSize', 14);

% Add a text label 'RNN'

% Make the axis lines and ticks visible and clean
box off, set(gca,'tickdir','out'),axis square
% Turn off the hold
hold off;