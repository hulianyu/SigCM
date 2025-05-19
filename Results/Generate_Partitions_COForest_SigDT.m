addpath([cd '/']);
addpath([cd '/Datasets']);
addpath([cd '/COForest']);
addpath([cd '/SigDT']);
%% Load Data sets
filename = char('lenses','lung-cancer','soybean-small','zoo','dna-promoter',...
    'hayes-roth','lymphography','heart-disease','solar-flare','primary-tumor',...
    'dermatology','house-votes','balance-scale','credit-approval','breast-cancer-wisconsin',...
    'mammographic-mass','tic-tac-toe','car');
rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf',...
    'Pt','De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce'};
TotalRun = 100;
%% Generate COForest Partitions
Partitions_COForest = cell(18,1);
for I=1:18
    disp(I);
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1);
    GT = X_data(:,1); %Ground Truth
    K = length(unique(GT)); %Cluster Number
    pi = zeros(N,TotalRun);
    parfor run=1:TotalRun
        pi(:,run) = COForest_main(X,K);
    end
    Partitions_COForest{I,1} = pi;
end
save('Partitions_COForest.mat', 'Partitions_COForest');

%% Generate SigDT Partitions
Partitions_SigDT = cell(18,1);
for I=1:18
    disp(I);
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    pi = SigDT_main(X);
    Partitions_SigDT{I,1} = pi;
end
save('Partitions_SigDT.mat', 'Partitions_SigDT');