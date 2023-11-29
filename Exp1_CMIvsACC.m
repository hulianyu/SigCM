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
        com_ACC_random(I,run) = ACC(GT, pi);
        com_CMI_random(I,run) = CMI(X,pi);
    end
end
com_ACC_random_mean = mean(com_ACC_random,2);
com_CMI_random_mean = mean(com_CMI_random,2);

%% Partitions: kmodes
load('Partitions_kmodes.mat');
RT = 100;
com_ACC_kmodes = zeros(18,RT);
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
        com_ACC_kmodes(I,run) = ACC(GT, pi);
        com_CMI_kmodes(I,run) = CMI(X,pi);
    end
end
com_ACC_kmodes_mean = mean(com_ACC_kmodes,2);
com_CMI_kmodes_mean = mean(com_CMI_kmodes,2);

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