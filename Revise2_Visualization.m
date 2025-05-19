addpath([cd '/']);
addpath([cd '/Datasets']);
addpath([cd '/InternalMetrics']);
Metrics = zeros(3,8);
%% Load Data sets
filename = char('zoo');
rowNames = {'Zo'};
%[4] Zo
load('Revise2_Visualization.mat')
%% Plot

figure; 

load('Partitions_random.mat');

% --------- 1. Partitions_random: Zo ---------
subplot(3,1,1);
Zo_random_labels = Partitions_random{4,1}(:,6);
%
all_metrics = ClusteringMeasure(GT, Zo_random_labels);
Metrics(1,1) = all_metrics(1);
Metrics(1,2) = all_metrics(2);
Metrics(1,3) = all_metrics(4);
Metrics(1,4) = all_metrics(7);
Metrics(1,5) = CMI(Zo_X,Zo_random_labels);
Metrics(1,6) = KMF(Zo_X,Zo_random_labels);
Metrics(1,7) = Entropy(Zo_X,Zo_random_labels);
Metrics(1,8) = CU(Zo_X,Zo_random_labels);
%
[~, ~, pi] = unique(Zo_random_labels);
maxK = max(pi);
remain = [];
for k =1:maxK
    clusterk = find(pi==k);
    pvals = SigCM_intra(Zo_X,clusterk);
    [~,~,list] = FWER_Control(pvals);
    if ~isempty(list)
        remain = [remain;list(:,1)];
    end
end
remain = unique(remain);

hold on
all_idx = 1:size(Zo_Y,1);
other_idx = setdiff(all_idx, remain);
if isempty(remain)
    scatter(Zo_Y(:,1), Zo_Y(:,2), 85, Zo_random_labels, 's', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.7);
else
    scatter(Zo_Y(remain,1), Zo_Y(remain,2), 60, Zo_random_labels(remain), 'o', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.7);
    scatter(Zo_Y(other_idx,1), Zo_Y(other_idx,2), 85, Zo_random_labels(other_idx), 's', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.7);
end
hold off
colormap(lines(numel(unique(Zo_random_labels))));
axis off;

load('Partitions_kmodes.mat')

% --------- 2. Partitions_kmodes: Zo ---------
subplot(3,1,2);
Zo_kmodes_labels = Partitions_kmodes{4,1}(:,6);
%
all_metrics = ClusteringMeasure(GT, Zo_kmodes_labels);
Metrics(2,1) = all_metrics(1);
Metrics(2,2) = all_metrics(2);
Metrics(2,3) = all_metrics(4);
Metrics(2,4) = all_metrics(7);
Metrics(2,5) = CMI(Zo_X,Zo_kmodes_labels);
Metrics(2,6) = KMF(Zo_X,Zo_kmodes_labels);
Metrics(2,7) = Entropy(Zo_X,Zo_kmodes_labels);
Metrics(2,8) = CU(Zo_X,Zo_kmodes_labels);
%
[~, ~, pi] = unique(Zo_kmodes_labels);
maxK = max(pi);
remain = [];
for k =1:maxK
    clusterk = find(pi==k);
    pvals = SigCM_intra(Zo_X,clusterk);
    [~,~,list] = FWER_Control(pvals);
    if ~isempty(list)
        remain = [remain;list(:,1)];
    end
end
remain = unique(remain);

hold on
all_idx = 1:size(Zo_Y,1);
other_idx = setdiff(all_idx, remain);
if isempty(remain)
    scatter(Zo_Y(:,1), Zo_Y(:,2), 85, Zo_kmodes_labels, 's', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.7);
else
    scatter(Zo_Y(remain,1), Zo_Y(remain,2), 60, Zo_kmodes_labels(remain), 'o', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.7);
    scatter(Zo_Y(other_idx,1), Zo_Y(other_idx,2), 85, Zo_kmodes_labels(other_idx), 's', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.7);
end
hold off
colormap(lines(numel(unique(Zo_kmodes_labels))));
axis off;

load('Partitions_CMS.mat')

% --------- 3. Partitions_CMS: Zo ---------
subplot(3,1,3);
Zo_CMS_labels = Partitions_CMS{4,1}(:,6);
%
all_metrics = ClusteringMeasure(GT, Zo_CMS_labels);
Metrics(3,1) = all_metrics(1);
Metrics(3,2) = all_metrics(2);
Metrics(3,3) = all_metrics(4);
Metrics(3,4) = all_metrics(7);
Metrics(3,5) = CMI(Zo_X,Zo_CMS_labels);
Metrics(3,6) = KMF(Zo_X,Zo_CMS_labels);
Metrics(3,7) = Entropy(Zo_X,Zo_CMS_labels);
Metrics(3,8) = CU(Zo_X,Zo_CMS_labels);
%
[~, ~, pi] = unique(Zo_CMS_labels);
maxK = max(pi);
remain = [];
for k =1:maxK
    clusterk = find(pi==k);
    pvals = SigCM_intra(Zo_X,clusterk);
    [~,~,list] = FWER_Control(pvals);
    if ~isempty(list)
        remain = [remain;list(:,1)];
    end
end
remain = unique(remain);

hold on
all_idx = 1:size(Zo_Y,1);
other_idx = setdiff(all_idx, remain);
if isempty(remain)
    scatter(Zo_Y(:,1), Zo_Y(:,2), 85, Zo_CMS_labels, 's', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.7);
else
    scatter(Zo_Y(remain,1), Zo_Y(remain,2), 60, Zo_CMS_labels(remain), 'o', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.7);
    scatter(Zo_Y(other_idx,1), Zo_Y(other_idx,2), 85, Zo_CMS_labels(other_idx), 's', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.7);
end
hold off
colormap(lines(numel(unique(Zo_CMS_labels))));
axis off;