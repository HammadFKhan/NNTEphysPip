%% Regression analysis of PC data with decodability
% Load in data
Spikes = makeSpikeGPFA(Spikes);
Spikes.GPFA.HitMiss.dat = [Spikes.GPFA.hit.dat,Spikes.GPFA.miss.dat];
for n = 1:length(IntanBehaviour.hitTrace)%+1:length(Spikes.GPFA.BaselineOpto.dat) %fix trials
    Spikes.GPFA.HitMiss.dat(n).trialId = n;
end
Spikes.GPFA.MIHitFA.dat = [Spikes.GPFA.MIHit.dat,Spikes.GPFA.MIFA.dat];
for n = length(IntanBehaviour.MIHitTrace)+1:length(Spikes.GPFA.MIHitFA.dat) %fix trials
    Spikes.GPFA.MIHitFA.dat(n).trialId = n;
end
%%%
addpath(genpath('C:\Users\khan332\Documents\GitHub\NeuralTraj'));
addpath(genpath('mat_results'));
if exist('mat_results','dir'),rmdir('mat_results','s'),end
[Spikes.GPFA.resultHit,Spikes.GPFA.seqTrainHit] = gpfaAnalysis(Spikes.GPFA.hit.dat,1); %Run index
[Spikes.GPFA.resultMiss,Spikes.GPFA.seqTrainMiss] = gpfaAnalysis(Spikes.GPFA.miss.dat,2); %Run index


%% Neural Trajectory Analysis
%IntanBehaviour.parameters = parameters;
%neuralTrajAnalysis(Spikes,Waves1,IntanBehaviour1);
[neuralDynamics,M1waveDynamics] = neuralTrajAnalysis2(Spikes,[],IntanBehaviour);

leverTrace = arrayfun(@(x) x.trace',IntanBehaviour.hitTrace,'uniformoutput',0);
leverTrace = vertcat(leverTrace{:});

reconstructedData = getPCALever(leverTrace);
leverTrace = reconstructedData;
leverTraceFull = reshape(leverTrace',[],1);
leverTraceFull = (leverTraceFull-min(leverTraceFull))./(max(leverTraceFull)-min(leverTraceFull));
leverTraceFull = leverTraceFull-min(leverTraceFull);
leverTraceFull = smoothdata(leverTraceFull,'gaussian',300);
leverTrace = (leverTrace-min(leverTrace,[],2))./(max(leverTrace,[],2)-min(leverTrace,[],2));
%%% Prime data up
dim = 15;
latentDynamicsSequence = squeeze(neuralDynamics.hitOnly.X(1:dim,:,:));
latentDynamicsSequence = reshape(latentDynamicsSequence,dim,[])';

% Compute mean and standard deviation for each PC
mu = mean(latentDynamicsSequence, 1); % Mean across rows (time × trials)
sigma = std(latentDynamicsSequence, 0, 1); % Std deviation across rows

% Standardize each column (PC)
standardizedlatentDynamicsSequence = ((latentDynamicsSequence - mu) ./ sigma)';
latentDynamicsSequence = reshape(standardizedlatentDynamicsSequence,dim,[]);

y = latentDynamicsSequence;
y_trial = reshape(latentDynamicsSequence,dim,250,[]);

% Check mean and variance for each PC
mean_check = mean(y, [2, 3]); % Mean across time and trials
var_check = var(y, [], [2, 3]); % Variance across time and trials

disp('Mean of standardized PCs:');
disp(mean(mean_check)); % Should be close to zero

disp('Variance of standardized PCs:');
disp(mean(var_check)); % Should be close to one


PC_avg = mean(y_trial,3);
downSampledLever = mean(leverTrace);
downSampledLever = resample(downSampledLever, size(PC_avg, 2), length(downSampledLever));
downSampledLeverFull = resample(leverTraceFull, size(y, 2), length(leverTraceFull));

%%% Build polynomial regression model

poly_mdl = fitlm(PC_avg', downSampledLever', 'quadratic');
disp(poly_mdl);
%figure,plot(poly_mdl)

% Fit polynomial regression

poly_mdlFull = fitlm(y', downSampledLeverFull', 'quadratic');
disp(poly_mdlFull);
figure,plot(poly_mdlFull)
%%% Extract predictors (PCs) and response from mdl
predictors = poly_mdlFull.Variables{:, 1:end-1}; % All columns except the last one (response)
response = poly_mdlFull.Variables{:, end};       % Last column is the response

%% Choose a PC to color by (e.g., first PC)
figure
color_by_pc = predictors(:, 1); % First principal component
scatter(response, poly_mdlFull.Fitted, [], color_by_pc, 'filled');
colorbar;
xlabel('Observed Lever Trace');
ylabel('Fitted Lever Trace');
title('Fitted vs Observed Lever Trace Colored by PC1');

figure
for n = 1:15
    color_by_pc = predictors(:, n); % First principal component
    subplot(4,4,n),scatter(response, poly_mdlFull.Fitted, [], color_by_pc, 'filled');
    colorbar;
    xlabel('Observed Lever Trace');
    ylabel('Fitted Lever Trace');
    title(['PC' num2str(n)]);
end
%%
residuals = poly_mdlFull.Residuals.Raw;

figure;
histogram(residuals, 'Normalization', 'probability', 'FaceColor', [0, 0.4470, 0.7410], 'EdgeColor', 'none');
xlabel('Residuals', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Probability', 'FontSize', 14, 'FontWeight', 'bold');
title('Residuals of Predicted Motor Behaviour', 'FontSize', 16, 'FontWeight', 'bold');

% Improve figure aesthetics for publication
set(gca, 'FontSize', 12); % Set axis font size
set(gca, 'LineWidth', 1); % Set axis line width
box on; % Add box around the plot
%% Compute R^2 for each individual PC
num_PCs = size(predictors, 2); % Number of PCs
R_squared_values = zeros(num_PCs, 1);

for i = 1:num_PCs
    mdl_single = fitlm(predictors(:, i), response,'quadratic'); % Fit model using one PC at a time
    R_squared_values(i) = mdl_single.Rsquared.Ordinary; % Extract R^2 value
end

figure,bar(R_squared_values)
%% Build regression decoder
% Example: 10-fold cross-validation
X = y';
Y = downSampledLeverFull;
cv = cvpartition(size(X, 1), 'KFold', 10); % X is your predictors (PCs)
mse_values = zeros(cv.NumTestSets, size(y,1));

for n = 1:size(y,1)
    X = y(n,:)';
    for i = 1:cv.NumTestSets
        % Training and test indices
        trainIdx = cv.training(i);
        testIdx = cv.test(i);

        % Training data
        X_train = X(trainIdx, :); % PCs for training
        y_train = Y(trainIdx);   % Lever trace for training

        % Test data
        X_test = X(testIdx, :);  % PCs for testing
        y_test = Y(testIdx);     % Lever trace for testing

        % Train the model (linear regression in this case)
        mdl = fitlm(X_train, y_train,"quadratic");

        % Predict on test data
        y_pred = predict(mdl, X_test);

        % Compute mean squared error for this fold on n dim
        mse_values(i,n) = mean((y_test - y_pred).^2);
    end
end
% Average MSE across folds
avg_mse = mean(mean(mse_values));
disp(['Cross-validated MSE: ', num2str(avg_mse)]);
%%% Compute R^2 for the fitted model
R_squared = mdl.Rsquared.Ordinary;
disp(['R^2: ', num2str(R_squared)]);
X = y';
% Compute correlation coefficient between predicted and observed lever traces
mdl = fitlm(X, Y,"quadratic");
correlation = corr(Y, predict(mdl, X));
disp(['Correlation Coefficient: ', num2str(correlation)]);
%%
% Plot observed vs predicted lever trace
figure;
scatter(Y, predict(mdl, X), 'filled');
xlabel('Observed Lever Trace');
ylabel('Predicted Lever Trace');
title('Decoder Performance: Observed vs Predicted Lever Trace');
%%
X = y';
% Predict lever trace using the trained decoder
predicted_lever_trace = predict(mdl, X); % X is your predictors (PCs)
% Time vector (adjust based on your data dimensions)
time =0:1/50:size(downSampledLeverFull,1)/50; % Assuming response corresponds to time points
time = time(2:end);
colors = [32/255 118/255 188/255; 233/255 14/255 139/255;24/255 158/255 118/255];

% Plot observed and predicted lever traces
figure;
subplot(211),plot(time, downSampledLeverFull, 'color',[0.5 0.5 0.5], 'LineWidth', 1); % Observed lever trace (blue)
hold on;
subplot(211),plot(time, predicted_lever_trace, 'r', 'LineWidth', 2); % Predicted lever trace (red dashed)
hold off;

% Add labels and legend
xlabel('Time');
ylabel('Lever Trace');
legend({'Observed Lever Trace', 'Predicted Lever Trace'});
title('Observed vs Predicted Lever Trace');
xlim([150 175])
subplot(2,1,2),plot(time,smoothdata(y(1,:),'gaussian',10),'color',colors(1,:),'linewidth',2),hold on
subplot(2,1,2),plot(time,smoothdata(y(2,:),'gaussian',10),'color',colors(2,:),'linewidth',2)
subplot(2,1,2),plot(time,smoothdata(y(3,:),'gaussian',10),'color',colors(3,:),'linewidth',2), box off
xlim([150 175])
xlabel('Time (s)')
ylabel('Latent Values (au)')

%%
% Compute error rate and accuracy over time
% we shift to correct for negative values
normalized_error = abs(predicted_lever_trace - downSampledLeverFull) / max(abs(downSampledLeverFull));
% Compute accuracy as 1 - normalized error
accuracy_trace = 1 - normalized_error; % Avoid division by zero
% Plot accuracy trace over time
figure;
subplot(211),plot(time, downSampledLeverFull, 'color',[0.5 0.5 0.5], 'LineWidth', 1); % Observed lever trace (blue)
hold on;
subplot(211),plot(time, predicted_lever_trace, 'r', 'LineWidth', 2); % Predicted lever trace (red dashed)
hold off;
subplot(212),plot(time, accuracy_trace, 'b', 'LineWidth', 2);hold on

xlabel('Time');
ylabel('Accuracy');
title('Accuracy of Predicted Lever Trace vs True Lever Trace Over Time');
grid on;
%%
% Compute lever speed (first derivative of lever trace)
lever_speed = diff(leverTraceFull) / 0.001; % dt is the time step between samples

% To match dimensions with other data, append a value (e.g., 0) at the start
lever_speed = [0; lever_speed];
downSampledLeverSpeed = resample(lever_speed, size(y, 2), length(lever_speed));
% Predictors (PCs): Time x PC dimensions
X = y'; % Replace with your actual PC data

% Responses (lever position and speed): Time x 2
Y = [downSampledLeverFull, downSampledLeverSpeed]; % Combine position and speed into one matrix

% Add intercept term to predictors
X_design = [ones(size(X, 1), 1), X]; % Add column of ones for intercept

% Fit multivariate regression model
[beta, sigma] = mvregress(X_design, Y);

% Display estimated coefficients
disp('Estimated Coefficients:');
disp(beta);

% Display covariance matrix of errors
disp('Error Covariance Matrix:');
disp(sigma);

% Predicted responses
Y_predicted = X_design * beta;

% Compute residuals
residuals = Y - Y_predicted;

% Plot residuals for each response variable
figure;
subplot(2, 1, 1);
plot(residuals(:, 1));
title('Residuals for Lever Position');
xlabel('Time');
ylabel('Residual');

subplot(2, 1, 2);
plot(residuals(:, 2));
title('Residuals for Lever Speed');
xlabel('Time');
ylabel('Residual');
%% Compute with new data
neuralDynamics = M1neuralDynamics(5).neuralDynamics;
IntanBehaviour = M1neuralDynamics(5).IntanBehaviour;
leverTrace = arrayfun(@(x) x.trace',IntanBehaviour.cueHitTrace,'uniformoutput',0);
leverTrace = vertcat(leverTrace{:});
reconstructedData = getPCALever(leverTrace);
leverTrace = reconstructedData;
leverTraceFull = reshape(leverTrace',[],1);
leverTraceFull = (leverTraceFull-min(leverTraceFull))./(max(leverTraceFull)-min(leverTraceFull));
leverTraceFull = leverTraceFull-min(leverTraceFull);
leverTraceFull = smoothdata(leverTraceFull,'gaussian',300);
%%% Prime data up
dim = 15;
latentDynamicsSequence = squeeze(neuralDynamics.hit.X(1:dim,:,:));
latentDynamicsSequence = reshape(latentDynamicsSequence,dim,[])';

% Compute mean and standard deviation for each PC
mu = mean(latentDynamicsSequence, 1); % Mean across rows (time × trials)
sigma = std(latentDynamicsSequence, 0, 1); % Std deviation across rows

% Standardize each column (PC)
standardizedlatentDynamicsSequence = ((latentDynamicsSequence - mu) ./ sigma)';
latentDynamicsSequence = reshape(standardizedlatentDynamicsSequence,dim,[]);

Xnew = latentDynamicsSequence';
Ynew = resample(leverTraceFull, size(latentDynamicsSequence, 2), length(leverTraceFull));
%%% Predict with new model
predicted_lever_tracenew = predict(mdl, Xnew); % X is your predictors (PCs)

% Compute MSE
mse = mean((predicted_lever_trace(1:length(predicted_lever_tracenew)) - predicted_lever_tracenew).^2);

% Compute correlation coefficient
correlation = corr(predicted_lever_trace(1:length(predicted_lever_tracenew)), predicted_lever_tracenew);

% Display results
disp(['MSE: ', num2str(mse)]);
disp(['Correlation Coefficient: ', num2str(correlation)]);
%%
time =0:1/50:size(Ynew,1)/50; % Assuming response corresponds to time points
time = time(2:end);

figure;
subplot(211),plot(time, Ynew, 'color',[0.5 0.5 0.5], 'LineWidth', 1); % Observed lever trace (blue)
hold on;
subplot(211),plot(time, predicted_lever_tracenew, 'r', 'LineWidth', 2); % Predicted lever trace (red dashed)
hold off;

% Add labels and legend
xlabel('Time');
ylabel('Lever Trace');
legend({'Observed Lever Trace', 'Predicted Lever Trace'});
title('Observed vs Predicted Lever Trace');
%%
win = 1:150;
figure
subplot(311),plot(squeeze(y_trial(1,win,:)),'color',[0 0 0 0.25])
subplot(312),plot(squeeze(y_trial(2,win,:)),'color',[0 0 0 0.25])
subplot(313),plot(squeeze(y_trial(3,win,:)),'color',[0 0 0 0.25])

%% LOCAL FUNCTIONS
function reconstructedData = getPCALever(data)
% Transpose the data
data = data';

% Perform PCA
[coeff, score, latent] = pca(data);

% Determine how many principal components to retain
cum_var = cumsum(latent ./ sum(latent));
n_components = find(cum_var >= 0.99, 1, 'first');
disp('Num of PCA components:')
disp(n_components)
% Select the first n_components principal components
selected_coeff = coeff(:, 1:n_components);
selected_score = score(:, 1:n_components);

% Calculate the mean of your original data
mean_data = mean(data);

% Reconstruct the data using the selected principal components
reconstructedData = selected_score * selected_coeff' + repmat(mean_data, size(data, 1), 1);
reconstructedData = reconstructedData';
end