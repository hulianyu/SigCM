addpath([cd '/']);
addpath([cd '/Datasets']);
%% Load Data sets
filename = char('lenses','lung-cancer','soybean-small','zoo','dna-promoter',...
    'hayes-roth','lymphography','heart-disease','solar-flare','primary-tumor',...
    'dermatology','house-votes','balance-scale','credit-approval','breast-cancer-wisconsin',...
    'mammographic-mass','tic-tac-toe','car');
rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf',...
    'Pt','De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce'};

%% Partitions: ADC
load('Partitions_ADC.mat');
RT = 100;
com_ACC_ADC = zeros(18,RT);
com_NMI_ADC = zeros(18,RT);
com_ARI_ADC = zeros(18,RT);
com_FSC_ADC = zeros(18,RT);
com_CMI_ADC = zeros(18,RT);
for I=1:18
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1); 
    GT = X_data(:,1); %Ground Truth
    %
    parfor run =1:RT
        disp('Partitions_ADC');
        disp(I);
        disp(run);
        pi = Partitions_ADC{I,1}(:,run);
        [~, ~, pi] = unique(pi);
        all_metrics = ClusteringMeasure(GT, pi);
        com_ACC_ADC(I,run) = all_metrics(1);
        com_NMI_ADC(I,run) = all_metrics(2);
        com_ARI_ADC(I,run) = all_metrics(4);
        com_FSC_ADC(I,run) = all_metrics(7);
        com_CMI_ADC(I,run) = CMI(X,pi);
    end
end
com_ACC_ADC_mean = mean(com_ACC_ADC,2);
com_NMI_ADC_mean = mean(com_NMI_ADC,2);
com_ARI_ADC_mean = mean(com_ARI_ADC,2);
com_FSC_ADC_mean = mean(com_FSC_ADC,2);
com_CMI_ADC_mean = mean(com_CMI_ADC,2);

%% Partitions: COForest
load('Partitions_COForest.mat');
RT = 100;
com_ACC_COForest = zeros(18,RT);
com_NMI_COForest = zeros(18,RT);
com_ARI_COForest = zeros(18,RT);
com_FSC_COForest = zeros(18,RT);
com_CMI_COForest = zeros(18,RT);
for I=1:18
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1); 
    GT = X_data(:,1); %Ground Truth
    %
    parfor run =1:RT
        disp('Partitions_COForest');
        disp(I);
        disp(run);
        pi = Partitions_COForest{I,1}(:,run);
        [~, ~, pi] = unique(pi);
        all_metrics = ClusteringMeasure(GT, pi);
        com_ACC_COForest(I,run) = all_metrics(1);
        com_NMI_COForest(I,run) = all_metrics(2);
        com_ARI_COForest(I,run) = all_metrics(4);
        com_FSC_COForest(I,run) = all_metrics(7);
        com_CMI_COForest(I,run) = CMI(X,pi);
    end
end
com_ACC_COForest_mean = mean(com_ACC_COForest,2);
com_NMI_COForest_mean = mean(com_NMI_COForest,2);
com_ARI_COForest_mean = mean(com_ARI_COForest,2);
com_FSC_COForest_mean = mean(com_FSC_COForest,2);
com_CMI_COForest_mean = mean(com_CMI_COForest,2);

%% Partitions: SigDT
load('Partitions_SigDT.mat');
com_ACC_SigDT_mean = zeros(18,1);
com_NMI_SigDT_mean = zeros(18,1);
com_ARI_SigDT_mean = zeros(18,1);
com_FSC_SigDT_mean = zeros(18,1);
com_CMI_SigDT_mean = zeros(18,1);
for I=1:18
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1);
    GT = X_data(:,1); %Ground Truth
    %
    disp('Partitions_SigDT');
    disp(I);
    pi = Partitions_SigDT{I,1};
    [~, ~, pi] = unique(pi);
    all_metrics = ClusteringMeasure(GT, pi);
    com_ACC_SigDT_mean(I,1) = all_metrics(1);
    com_NMI_SigDT_mean(I,1) = all_metrics(2);
    com_ARI_SigDT_mean(I,1) = all_metrics(4);
    com_FSC_SigDT_mean(I,1) = all_metrics(7);
    com_CMI_SigDT_mean(I,1) = CMI(X,pi);
end

%% save 'com_' vars 
vars = who('com_*');
save('Revise2_Exp1_CMIvsACC_NMI_ARI_Fscore.mat', vars{:});