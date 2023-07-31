clear
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.mat'));
L = length(directory);
for i = 2
    disp(['Parsing: ' num2str(directory(i).name)])
    S = load(fullfile(directory(i).folder,directory(i).name));
    StrName=fieldnames(S);
    StrName = StrName{1};
%     [P_c{i},I_c{i},H_c(i,1),H_e(i,1),sizeE{i},sizeEdge{i}] = EntropyParser(S.(StrName));
    [theta{i},beta{i},gamma{i}] = powerSpecStat(S.(StrName));
end 

thetaCon2 = vertcat(theta{:});betaCon2 = vertcat(beta{:});gammaCon2 = vertcat(gamma{:});
b1 = mean(betaCon1);
betaM = cellfun(@mean,allBeta);
%%
for i = 2:6
    data = allBeta{i}-bT;
    figure,
    boxplot(data-bT,i,'PlotStyle','compact'),ylim([-40 40]),box off
end