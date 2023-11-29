addpath([cd '/']);
addpath([cd '/Datasets']);
%% Load Data sets
filename = char('lenses','lung-cancer','soybean-small','zoo','dna-promoter',...
    'hayes-roth','lymphography','heart-disease','solar-flare','primary-tumor',...
    'dermatology','house-votes','balance-scale','credit-approval','breast-cancer-wisconsin',...
    'mammographic-mass','tic-tac-toe','car');
rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf',...
    'Pt','De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce'};
TotalRun = 100;
%% Generate Random Partitions
Partitions_random = cell(18,1);
for I=1:18
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    N = size(X_data,1);
    GT = X_data(:,1); %Ground Truth
    pi = zeros(N,TotalRun);
    for run=1:TotalRun
        pi(:,run) = Pi_randperm(GT);
    end
    Partitions_random{I,1} = pi;
end
% save('Partitions_random.mat', 'Partitions_random');


%% Generate k-modes Partitions
Partitions_kmodes = cell(18,1);
for I=1:18
    disp(I);
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    N = size(X,1);
    GT = X_data(:,1); %Ground Truth
    K = length(unique(GT)); %Cluster Number
    pi = zeros(N,TotalRun);
    parfor run=1:TotalRun
        pi(:,run) = kmode_random(X,K);
    end
    Partitions_kmodes{I,1} = pi;
end
% save('Partitions_kmodes.mat', 'Partitions_kmodes');

%%%%%%%%%%%%%%%%%%%%%%% Functions for generation %%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function pi = Pi_randperm(GT)
% Generate random partition pi based on ground-truth
% randperm(N)
N = length(GT);
pi = GT(randperm(N),1); %randomly permuting the ground-truth vector
end
%%
function obs_cluster = kmode_random(data, k)
m = size(data,2);
n = size(data,1);
data(:,end+1) = 0;
% Randomly select k observations from the data to initialize the modes
randIndices = randperm(n, k);
mode = data(randIndices,1:m);
% mode = data(1:k,1:m) ;
while 1
    indicator = 1;
    for i = 1:n
        d = zeros(1,k);
        for j = 1:k
            d(1,j) = dist_cate(data(i,1:m), mode(j,1:m));
        end
        [~, ind] = min(d);
        indicator = indicator & (ind == data(i,end));
        data(i, end) = ind;
        cluster = data(data(:,end) == ind, :);
        mode(ind,:) = get_mode0(cluster(:,1:m));
    end
    if indicator == 1
        break;
    end
end
obs_cluster = data(:,end);
end
%%
function mode = get_mode0(data)
% the purpose is to obtain mode(s) of the data set, by choosing most frequent values of each column
m = length(data(1,:));
mode = zeros(1,m);
for i = 1: m
    sort_col = sort(data(:,i))';
    run_length = diff([ 0 find(sort_col(1:end-1) ~= sort_col(2:end)) length(sort_col) ]);
    uni_col = unique(sort_col);
    [~, ind] = max(run_length);
    %ind = find(run_length == val);
    mode(i) = uni_col(ind);
end
end
%%
function distance = dist_cate (a,b)
% purpose: calculate distance between two vector with distance defined by summation of
%          distance of each attribute, where distance of attribute is the indicator
%          function of non-equality.
% record of revisions:
%     date               programmer              description of change
% -----------        -----------------          ------------------------
% June 10,2003           Peng Zhang                 Original code
%
% define variables:
% a           -- vector a
% b           -- vector b
% distance    -- returned distance between a and b
%
% check whether a and b have the same length, if not set distance -1
if length(a) ~= length(b)
    distance = -1;
else
    distance = 0;
    for i = 1:length(a)
        distance = distance + (a(i) ~= b(i));
    end
end
end