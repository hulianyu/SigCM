addpath([cd '/']);
addpath([cd '/Datasets']);
%% Load Data sets
filename = char('lenses','lung-cancer','soybean-small','zoo','dna-promoter',...
    'hayes-roth','lymphography','heart-disease','solar-flare','primary-tumor',...
    'dermatology','house-votes','balance-scale','credit-approval','breast-cancer-wisconsin',...
    'mammographic-mass','tic-tac-toe','car');
rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf',...
    'Pt','De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce'};
%% Partitions: random
load('Partitions_random.mat');
RT = 100;
com_ACC_random = zeros(18,RT);
com_NMI_random = zeros(18,RT);
com_ARI_random = zeros(18,RT);
com_FSC_random = zeros(18,RT);
com_CMI_random = zeros(18,RT);
for I=1:18
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1); 
    GT = X_data(:,1); %Ground Truth
    %
    parfor run =1:RT
        disp('Partitions_random');
        disp(I);
        disp(run);
        pi = Partitions_random{I,1}(:,run);
        [~, ~, pi] = unique(pi);
        all_metrics = ClusteringMeasure(GT, pi);
        com_ACC_random(I,run) = all_metrics(1);
        com_NMI_random(I,run) = all_metrics(2);
        com_ARI_random(I,run) = all_metrics(4);
        com_FSC_random(I,run) = all_metrics(7);
        com_CMI_random(I,run) = CMI(X,pi);
    end
end
com_ACC_random_mean = mean(com_ACC_random,2);
com_NMI_random_mean = mean(com_NMI_random,2);
com_ARI_random_mean = mean(com_ARI_random,2);
com_FSC_random_mean = mean(com_FSC_random,2);
com_CMI_random_mean = mean(com_CMI_random,2);

%% Partitions: kmodes
load('Partitions_kmodes.mat');
RT = 100;
com_ACC_kmodes = zeros(18,RT);
com_NMI_kmodes = zeros(18,RT);
com_ARI_kmodes = zeros(18,RT);
com_FSC_kmodes = zeros(18,RT);
com_CMI_kmodes = zeros(18,RT);
for I=1:18
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1); 
    GT = X_data(:,1); %Ground Truth
    %
    parfor run =1:RT
        disp('Partitions_kmodes');
        disp(I);
        disp(run);
        pi = Partitions_kmodes{I,1}(:,run);
        [~, ~, pi] = unique(pi);
        all_metrics = ClusteringMeasure(GT, pi);
        com_ACC_kmodes(I,run) = all_metrics(1);
        com_NMI_kmodes(I,run) = all_metrics(2);
        com_ARI_kmodes(I,run) = all_metrics(4);
        com_FSC_kmodes(I,run) = all_metrics(7);
        com_CMI_kmodes(I,run) = CMI(X,pi);
    end
end
com_ACC_kmodes_mean = mean(com_ACC_kmodes,2);
com_NMI_kmodes_mean = mean(com_NMI_kmodes,2);
com_ARI_kmodes_mean = mean(com_ARI_kmodes,2);
com_FSC_kmodes_mean = mean(com_FSC_kmodes,2);
com_CMI_kmodes_mean = mean(com_CMI_kmodes,2);

%% Partitions: CDCDR
load('Partitions_CDCDR.mat');
RT = 100;
com_ACC_CDCDR = zeros(18,RT);
com_NMI_CDCDR = zeros(18,RT);
com_ARI_CDCDR = zeros(18,RT);
com_FSC_CDCDR = zeros(18,RT);
com_CMI_CDCDR = zeros(18,RT);
for I=1:18
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1); 
    GT = X_data(:,1); %Ground Truth
    %
    parfor run =1:RT
        disp('Partitions_CDCDR');
        disp(I);
        disp(run);
        pi = Partitions_CDCDR{I,1}(:,run);
        [~, ~, pi] = unique(pi);
        all_metrics = ClusteringMeasure(GT, pi);
        com_ACC_CDCDR(I,run) = all_metrics(1);
        com_NMI_CDCDR(I,run) = all_metrics(2);
        com_ARI_CDCDR(I,run) = all_metrics(4);
        com_FSC_CDCDR(I,run) = all_metrics(7);
        com_CMI_CDCDR(I,run) = CMI(X,pi);
    end
end
com_ACC_CDCDR_mean = mean(com_ACC_CDCDR,2);
com_NMI_CDCDR_mean = mean(com_NMI_CDCDR,2);
com_ARI_CDCDR_mean = mean(com_ARI_CDCDR,2);
com_FSC_CDCDR_mean = mean(com_FSC_CDCDR,2);
com_CMI_CDCDR_mean = mean(com_CMI_CDCDR,2);

%% Partitions: CMS
load('Partitions_CMS.mat');
RT = 100;
com_ACC_CMS = zeros(18,RT);
com_NMI_CMS = zeros(18,RT);
com_ARI_CMS = zeros(18,RT);
com_FSC_CMS = zeros(18,RT);
com_CMI_CMS = zeros(18,RT);
for I=1:18
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1); 
    GT = X_data(:,1); %Ground Truth
    %
    parfor run =1:RT
        disp('Partitions_CMS');
        disp(I);
        disp(run);
        pi = Partitions_CMS{I,1}(:,run);
        [~, ~, pi] = unique(pi);
        all_metrics = ClusteringMeasure(GT, pi);
        com_ACC_CMS(I,run) = all_metrics(1);
        com_NMI_CMS(I,run) = all_metrics(2);
        com_ARI_CMS(I,run) = all_metrics(4);
        com_FSC_CMS(I,run) = all_metrics(7);
        com_CMI_CMS(I,run) = CMI(X,pi);
    end
end
com_ACC_CMS_mean = mean(com_ACC_CMS,2);
com_NMI_CMS_mean = mean(com_NMI_CMS,2);
com_ARI_CMS_mean = mean(com_ARI_CMS,2);
com_FSC_CMS_mean = mean(com_FSC_CMS,2);
com_CMI_CMS_mean = mean(com_CMI_CMS,2);

%% Partitions: Het2Hom
load('Partitions_Het2Hom.mat');
RT = 100;
com_ACC_Het2Hom = zeros(18,RT);
com_NMI_Het2Hom = zeros(18,RT);
com_ARI_Het2Hom = zeros(18,RT);
com_FSC_Het2Hom = zeros(18,RT);
com_CMI_Het2Hom = zeros(18,RT);
for I=1:18
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1); 
    GT = X_data(:,1); %Ground Truth
    %
    parfor run =1:RT
        disp('Partitions_Het2Hom');
        disp(I);
        disp(run);
        pi = Partitions_Het2Hom{I,1}(:,run);
        [~, ~, pi] = unique(pi);
        all_metrics = ClusteringMeasure(GT, pi);
        com_ACC_Het2Hom(I,run) = all_metrics(1);
        com_NMI_Het2Hom(I,run) = all_metrics(2);
        com_ARI_Het2Hom(I,run) = all_metrics(4);
        com_FSC_Het2Hom(I,run) = all_metrics(7);
        com_CMI_Het2Hom(I,run) = CMI(X,pi);
    end
end
com_ACC_Het2Hom_mean = mean(com_ACC_Het2Hom,2);
com_NMI_Het2Hom_mean = mean(com_NMI_Het2Hom,2);
com_ARI_Het2Hom_mean = mean(com_ARI_Het2Hom,2);
com_FSC_Het2Hom_mean = mean(com_FSC_Het2Hom,2);
com_CMI_Het2Hom_mean = mean(com_CMI_Het2Hom,2);

%% Partitions: ground-truth
RT = 100;
com_CMI_GT = zeros(18,1);
for I=1:18
    disp('ground-truth');
    disp(I);
    %% Choose I-th ODS
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1); 
    GT = X_data(:,1); %Ground Truth
    [~, ~, GT] = unique(GT);
    com_CMI_GT(I,1) = CMI(X,GT);
end

%% save 'com_' vars 
vars = who('com_*');
% save('Exp1_CMIvsACC_NMI_ARI_Fscore.mat', vars{:});