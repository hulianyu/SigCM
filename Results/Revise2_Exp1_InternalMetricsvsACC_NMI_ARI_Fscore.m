addpath([cd '/']);
addpath([cd '/Datasets']);
addpath([cd '/InternalMetrics']);
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
com_KMF_random = zeros(18,RT);
com_Entropy_random = zeros(18,RT);
com_CU_random = zeros(18,RT);
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
        com_KMF_random(I,run) = KMF(X,pi);
        com_Entropy_random(I,run) = Entropy(X,pi);
        com_CU_random(I,run) = CU(X,pi);
    end
end
com_ACC_random_mean = mean(com_ACC_random,2);
com_NMI_random_mean = mean(com_NMI_random,2);
com_ARI_random_mean = mean(com_ARI_random,2);
com_FSC_random_mean = mean(com_FSC_random,2);
com_KMF_random_mean = mean(com_KMF_random,2);
com_Entropy_random_mean = mean(com_Entropy_random,2);
com_CU_random_mean = mean(com_CU_random,2);
%% Partitions: kmodes
load('Partitions_kmodes.mat');
RT = 100;
com_ACC_kmodes = zeros(18,RT);
com_NMI_kmodes = zeros(18,RT);
com_ARI_kmodes = zeros(18,RT);
com_FSC_kmodes = zeros(18,RT);
com_KMF_kmodes = zeros(18,RT);
com_Entropy_kmodes = zeros(18,RT);
com_CU_kmodes = zeros(18,RT);
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
        com_KMF_kmodes(I,run) = KMF(X,pi);
        com_Entropy_kmodes(I,run) = Entropy(X,pi);
        com_CU_kmodes(I,run) = CU(X,pi);
    end
end
com_ACC_kmodes_mean = mean(com_ACC_kmodes,2);
com_NMI_kmodes_mean = mean(com_NMI_kmodes,2);
com_ARI_kmodes_mean = mean(com_ARI_kmodes,2);
com_FSC_kmodes_mean = mean(com_FSC_kmodes,2);
com_KMF_kmodes_mean = mean(com_KMF_kmodes,2);
com_Entropy_kmodes_mean = mean(com_Entropy_kmodes,2);
com_CU_kmodes_mean = mean(com_CU_kmodes,2);

%% Partitions: CDCDR
load('Partitions_CDCDR.mat');
RT = 100;
com_ACC_CDCDR = zeros(18,RT);
com_NMI_CDCDR = zeros(18,RT);
com_ARI_CDCDR = zeros(18,RT);
com_FSC_CDCDR = zeros(18,RT);
com_KMF_CDCDR = zeros(18,RT);
com_Entropy_CDCDR = zeros(18,RT);
com_CU_CDCDR = zeros(18,RT);
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
        com_KMF_CDCDR(I,run) = KMF(X,pi);
        com_Entropy_CDCDR(I,run) = Entropy(X,pi);
        com_CU_CDCDR(I,run) = CU(X,pi);
    end
end
com_ACC_CDCDR_mean = mean(com_ACC_CDCDR,2);
com_NMI_CDCDR_mean = mean(com_NMI_CDCDR,2);
com_ARI_CDCDR_mean = mean(com_ARI_CDCDR,2);
com_FSC_CDCDR_mean = mean(com_FSC_CDCDR,2);
com_KMF_CDCDR_mean = mean(com_KMF_CDCDR,2);
com_Entropy_CDCDR_mean = mean(com_Entropy_CDCDR,2);
com_CU_CDCDR_mean = mean(com_CU_CDCDR,2);

%% Partitions: CMS
load('Partitions_CMS.mat');
RT = 100;
com_ACC_CMS = zeros(18,RT);
com_NMI_CMS = zeros(18,RT);
com_ARI_CMS = zeros(18,RT);
com_FSC_CMS = zeros(18,RT);
com_KMF_CMS = zeros(18,RT);
com_Entropy_CMS = zeros(18,RT);
com_CU_CMS = zeros(18,RT);
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
        com_KMF_CMS(I,run) = KMF(X,pi);
        com_Entropy_CMS(I,run) = Entropy(X,pi);
        com_CU_CMS(I,run) = CU(X,pi);
    end
end
com_ACC_CMS_mean = mean(com_ACC_CMS,2);
com_NMI_CMS_mean = mean(com_NMI_CMS,2);
com_ARI_CMS_mean = mean(com_ARI_CMS,2);
com_FSC_CMS_mean = mean(com_FSC_CMS,2);
com_KMF_CMS_mean = mean(com_KMF_CMS,2);
com_Entropy_CMS_mean = mean(com_Entropy_CMS,2);
com_CU_CMS_mean = mean(com_CU_CMS,2);

%% Partitions: Het2Hom
load('Partitions_Het2Hom.mat');
RT = 100;
com_ACC_Het2Hom = zeros(18,RT);
com_NMI_Het2Hom = zeros(18,RT);
com_ARI_Het2Hom = zeros(18,RT);
com_FSC_Het2Hom = zeros(18,RT);
com_KMF_Het2Hom = zeros(18,RT);
com_Entropy_Het2Hom = zeros(18,RT);
com_CU_Het2Hom = zeros(18,RT);
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
        com_KMF_Het2Hom(I,run) = KMF(X,pi);
        com_Entropy_Het2Hom(I,run) = Entropy(X,pi);
        com_CU_Het2Hom(I,run) = CU(X,pi);
    end
end
com_ACC_Het2Hom_mean = mean(com_ACC_Het2Hom,2);
com_NMI_Het2Hom_mean = mean(com_NMI_Het2Hom,2);
com_ARI_Het2Hom_mean = mean(com_ARI_Het2Hom,2);
com_FSC_Het2Hom_mean = mean(com_FSC_Het2Hom,2);
com_KMF_Het2Hom_mean = mean(com_KMF_Het2Hom,2);
com_Entropy_Het2Hom_mean = mean(com_Entropy_Het2Hom,2);
com_CU_Het2Hom_mean = mean(com_CU_Het2Hom,2);

%% Partitions: ADC
load('Partitions_ADC.mat');
RT = 100;
com_ACC_ADC = zeros(18,RT);
com_NMI_ADC = zeros(18,RT);
com_ARI_ADC = zeros(18,RT);
com_FSC_ADC = zeros(18,RT);
com_KMF_ADC = zeros(18,RT);
com_Entropy_ADC = zeros(18,RT);
com_CU_ADC = zeros(18,RT);
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
        com_KMF_ADC(I,run) = KMF(X,pi);
        com_Entropy_ADC(I,run) = Entropy(X,pi);
        com_CU_ADC(I,run) = CU(X,pi);
    end
end
com_ACC_ADC_mean = mean(com_ACC_ADC,2);
com_NMI_ADC_mean = mean(com_NMI_ADC,2);
com_ARI_ADC_mean = mean(com_ARI_ADC,2);
com_FSC_ADC_mean = mean(com_FSC_ADC,2);
com_KMF_ADC_mean = mean(com_KMF_ADC,2);
com_Entropy_ADC_mean = mean(com_Entropy_ADC,2);
com_CU_ADC_mean = mean(com_CU_ADC,2);

%% Partitions: COForest
load('Partitions_COForest.mat');
RT = 100;
com_ACC_COForest = zeros(18,RT);
com_NMI_COForest = zeros(18,RT);
com_ARI_COForest = zeros(18,RT);
com_FSC_COForest = zeros(18,RT);
com_KMF_COForest = zeros(18,RT);
com_Entropy_COForest = zeros(18,RT);
com_CU_COForest = zeros(18,RT);
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
        com_KMF_COForest(I,run) = KMF(X,pi);
        com_Entropy_COForest(I,run) = Entropy(X,pi);
        com_CU_COForest(I,run) = CU(X,pi);
    end
end
com_ACC_COForest_mean = mean(com_ACC_COForest,2);
com_NMI_COForest_mean = mean(com_NMI_COForest,2);
com_ARI_COForest_mean = mean(com_ARI_COForest,2);
com_FSC_COForest_mean = mean(com_FSC_COForest,2);
com_KMF_COForest_mean = mean(com_KMF_COForest,2);
com_Entropy_COForest_mean = mean(com_Entropy_COForest,2);
com_CU_COForest_mean = mean(com_CU_COForest,2);

%% Partitions: SigDT
load('Partitions_SigDT.mat');
com_ACC_SigDT_mean = zeros(18,1);
com_NMI_SigDT_mean = zeros(18,1);
com_ARI_SigDT_mean = zeros(18,1);
com_FSC_SigDT_mean = zeros(18,1);
com_KMF_SigDT_mean = zeros(18,1);
com_Entropy_SigDT_mean = zeros(18,1);
com_CU_SigDT_mean = zeros(18,1);
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
    com_KMF_SigDT_mean(I,1) = KMF(X,pi);
    com_Entropy_SigDT_mean(I,1) = Entropy(X,pi);
    com_CU_SigDT_mean(I,1) = CU(X,pi);
end

%% save 'com_' vars 
vars = who('com_*');
%save('Revise2_Exp1_InternalMetricsvsACC_NMI_ARI_Fscore.mat', vars{:});