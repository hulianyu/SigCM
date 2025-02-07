addpath([cd '/']);
addpath([cd '/Datasets']);
addpath([cd '/CDC_DR']);
addpath([cd '/CMS']);
addpath([cd '/CMS/Ncut']);
addpath([cd '/Het2Hom']);
%% Load Data sets
filename = char('lenses','lung-cancer','soybean-small','zoo','dna-promoter',...
    'hayes-roth','lymphography','heart-disease','solar-flare','primary-tumor',...
    'dermatology','house-votes','balance-scale','credit-approval','breast-cancer-wisconsin',...
    'mammographic-mass','tic-tac-toe','car');
rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf',...
    'Pt','De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce'};
TotalRun = 100;
%% Generate CDCDR Partitions
Partitions_CDCDR = cell(18,1);
for I=1:18
    disp(I);
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1);
    GT = X_data(:,1); %Ground Truth
    K = length(unique(GT)); %Cluster Number
    pi = zeros(N,TotalRun);
    parfor run=1:TotalRun
        pi(:,run) = CDC_DR_AE(X,K);
    end
    Partitions_CDCDR{I,1} = pi;
end
% save('Partitions_CDCDR.mat', 'Partitions_CDCDR');

%% Generate CMS Partitions
Partitions_CMS = cell(18,1);
for I=1:18
    disp(I);
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1);
    GT = X_data(:,1); %Ground Truth
    K = length(unique(GT)); %Cluster Number
    pi = zeros(N,TotalRun);
    parfor run=1:TotalRun
        pi(:,run) = CMS(X,K);
    end
    Partitions_CMS{I,1} = pi;
end
% save('Partitions_CMS.mat', 'Partitions_CMS');


%% Generate Het2Hom Partitions
Partitions_Het2Hom = cell(18,1);
for I=1:18
    disp(I);
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    [N,M] = size(X);
    ordinal_num = 0;
    nominal_num = M-ordinal_num;
    GT = X_data(:,1); %Ground Truth
    K = length(unique(GT)); %Cluster Number
    pi = zeros(N,TotalRun);
    parfor run=1:TotalRun
        pi(:,run) = Het2Hom_Clustering(X,K,nominal_num,ordinal_num)
    end
    Partitions_Het2Hom{I,1} = pi;
end
% save('Partitions_Het2Hom.mat', 'Partitions_Het2Hom');