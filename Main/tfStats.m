function stats = tfStats(TimeFreq)
L23 = TimeFreq.tf.depth.L23;
L4 = TimeFreq.tf.depth.L4;
L5 = TimeFreq.tf.depth.L5;

layersTheta = [L23.theta.itpc;L4.theta.itpc;L5.theta.itpc];
xTheta = [repmat({'L2/3 theta'},length(L23.theta.itpc), 1); repmat({'L4-theta'},length(L4.theta.itpc), 1); repmat({'L5-theta'},length(L5.theta.itpc), 1)];
layersBeta = [L23.beta.itpc;L4.beta.itpc;L5.beta.itpc];
xBeta = [repmat({'L2/3-beta'},length(L23.beta.itpc), 1); repmat({'L4-beta'},length(L4.beta.itpc), 1); repmat({'L5-beta'},length(L5.beta.itpc), 1)];
layersGamma = [L23.gamma.itpc;L4.gamma.itpc;L5.gamma.itpc];
xGamma = [repmat({'L2/3-gamma'},length(L23.gamma.itpc), 1); repmat({'L4-gamma'},length(L4.gamma.itpc), 1); repmat({'L5-gamma'},length(L5.gamma.itpc), 1)];

layers = [layersTheta; layersBeta; layersGamma];
x = [xTheta; xBeta; xGamma];

figure,boxplot(layers,x,'PlotStyle','compact'),hold on

x = [1.2*ones(length(L23.theta.itpc),1);2.2*ones(length(L4.theta.itpc),1);3.2*ones(length(L5.theta.itpc),1);...
    4.2*ones(length(L23.beta.itpc),1);5.2*ones(length(L4.beta.itpc),1);6.2*ones(length(L5.beta.itpc),1);...
    7.2*ones(length(L23.gamma.itpc),1);8.2*ones(length(L4.gamma.itpc),1);9.2*ones(length(L5.gamma.itpc),1)];
scatter(x,layers,'filled','k')


layers23 = [L23.theta.itpc;L23.beta.itpc;L23.gamma.itpc];
x23 = [repmat({'L2/3-theta'},length(L23.theta.itpc), 1); repmat({'L23-beta'},length(L23.beta.itpc), 1); repmat({'L23-gamma'},length(L23.gamma.itpc), 1)];
layers4 = [L4.theta.itpc;L4.beta.itpc;L4.gamma.itpc];
x4 = [repmat({'L4-theta'},length(L4.theta.itpc), 1); repmat({'L4-beta'},length(L4.beta.itpc), 1); repmat({'L4-gamma'},length(L4.gamma.itpc), 1)];
layers5 = [L5.theta.itpc;L5.beta.itpc;L5.gamma.itpc];
x5 = [repmat({'L5-theta'},length(L5.theta.itpc), 1); repmat({'L5-beta'},length(L5.beta.itpc), 1); repmat({'L5-gamma'},length(L5.gamma.itpc), 1)];

freq = [layers23; layers4; layers5];
xl = [x23; x4; x5];

figure,boxplot(freq,xl,'PlotStyle','compact')
% Export stat structure
stats.layersTheta = layersTheta;
stats.xTheta = xTheta;
stats.xBeta = xBeta;
stats.xGamma = xGamma;
stats.layers = layers;
stats.x = x;

stats.layers23 = layers23;
stats.layers4 = layers4;
stats.layers5 = layers5;
stats.freq = freq;
stats.xl = xl;
